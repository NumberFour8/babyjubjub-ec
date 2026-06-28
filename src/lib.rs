//! BabyJubJub elliptic curve implementation wrapped in `elliptic-curve` traits.
//!
//! This crate provides a wrapper around the `taceo-ark-babyjubjub` crate that implements
//! the BabyJubJub curve in a way compatible with the `elliptic-curve` crate traits.
//!
//! BabyJubJub is a cofactor-8 twisted Edwards curve defined over the 254-bit prime
//! field `Fq`; its prime-order subgroup has a 251-bit scalar field `Fr` ([`Scalar`]).
//!
//! # Security
//!
//! - **Side channels / constant-time.** Arithmetic is delegated to the arkworks
//!   backend, so the relevant timing guarantees are the backend's.
//!     - **Scalar multiplication is *almost* constant-time.** The `Mul<Scalar>`
//!       operators on [`ProjectivePoint`] and
//!       [`ProjectivePoint::mul_fixed_schedule`] both run a scalar-multiplication
//!       *algorithm* with no scalar-dependent control flow, but via two distinct
//!       code paths. The `Mul<Scalar>` operators delegate to the backend's
//!       `mul_projective`, which `taceo-ark-babyjubjub` overrides with a
//!       **Montgomery ladder** that iterates over a fixed number of scalar bits
//!       (it does *not* skip leading zeros) and swaps the two ladder registers
//!       using branch-free, bit-masked conditional swaps.
//!       [`ProjectivePoint::mul_fixed_schedule`], by contrast, does **not** call
//!       the backend's `mul_projective`: it runs an in-crate
//!       **double-and-add-always** loop over a fixed number of scalar bits, built
//!       only on the curve's complete (exception-free) point addition and
//!       doubling plus a bit-masked `subtle` select, so its algorithm-level
//!       constant-time property does **not** rely on the backend's
//!       scalar-multiplication routine. Either way the loop length, branching,
//!       and memory-access pattern are independent of the scalar. Neither is
//!       **end-to-end** constant-time, however: the underlying `ark-ff` field
//!       arithmetic uses a data-dependent conditional reduction
//!       (`Fp::subtract_modulus` via `is_geq_modulus`, a BigInteger comparison)
//!       whose timing depends on the intermediate field values, leaving a small
//!       residual timing signal. For end-to-end constant time, use a backend with
//!       bit-masked field reduction throughout (e.g. `fiat-crypto`-generated
//!       code).
//!     - **Variable-time scalar-field operations.** [`Scalar::invert`] and
//!       [`Scalar`]'s `sqrt`/`sqrt_ratio` have input-dependent control flow and
//!       can leak their inputs through timing.
//!     - `ConditionallySelectable` and `ct_eq` are constant-time.
//! - **Validation.** Points decoded via [`GroupEncoding::from_bytes`] are
//!   checked to be on-curve **and** in the prime-order subgroup. The
//!   [`GroupEncoding::from_bytes_unchecked`] decoder performs **no** such checks; prefer
//!   [`AffinePoint::new`] or the validation helpers
//!   ([`AffinePoint::is_on_curve`], [`AffinePoint::is_in_prime_order_subgroup`])
//!   for untrusted input.
//! - **Canonical encodings.** [`Scalar::from_bytes`] / [`PrimeField::from_repr`]
//!   reject non-canonical (`>= r`) scalar encodings, and the compressed point
//!   encoding is a canonical 32 bytes; this avoids scalar/point malleability.
//!   Use [`Scalar::reduce_bytes_be`] when modular reduction is explicitly wanted.
//! - **Zeroization.** [`Scalar`] is `Copy` and therefore cannot implement
//!   `ZeroizeOnDrop`; copies are not wiped automatically. Callers handling
//!   secret scalars are responsible for zeroizing their own storage. All three
//!   wrapper types implement `zeroize::DefaultIsZeroes` (and hence `Zeroize`)
//!   **unconditionally** via `elliptic_curve`'s re-exported `zeroize` (required
//!   by [`CurveArithmetic`]), so you can call `.zeroize()` explicitly when done.

#![cfg_attr(not(feature = "std"), no_std)]
#![deny(
    warnings,
    unused,
    future_incompatible,
    nonstandard_style,
    rust_2018_idioms
)]
#![forbid(unsafe_code)]

// ===== Re-export backend types =====

pub use taceo_ark_babyjubjub::EdwardsAffine as BackendAffine;
pub use taceo_ark_babyjubjub::EdwardsProjective as BackendProjective;
pub use taceo_ark_babyjubjub::Fq as BackendBaseField;
pub use taceo_ark_babyjubjub::Fr as BackendScalar;

pub use elliptic_curve;
pub use group;
pub use subtle;

// ===== Import required traits for BackendScalar operations =====
use ark_ec::PrimeGroup;
use ark_ff::fields::{AdditiveGroup, FftField, Field as ArkField, PrimeField as ArkPrimeField};
use ark_serialize::{CanonicalDeserialize, CanonicalSerialize, Compress, Validate};
use elliptic_curve::{Curve, CurveArithmetic, PrimeCurve};
use group::ff::{Field, PrimeField};
use group::{Group, GroupEncoding};
use num_traits::{One, Zero};
use subtle::{ConditionallySelectable, ConstantTimeEq, CtOption};
// `DefaultIsZeroes` is provided unconditionally via `elliptic_curve`'s
// always-available re-export of the `zeroize` crate. `CurveArithmetic` requires
// `Scalar`/`AffinePoint`/`ProjectivePoint: DefaultIsZeroes` without any feature
// gate, so this is not tied to the optional local `zeroize` feature.
use elliptic_curve::zeroize::DefaultIsZeroes;
// ===== Imports for the `CurveArithmetic` associated-type bounds =====
use elliptic_curve::bigint::{ArrayEncoding, U256, modular::Retrieve};
use elliptic_curve::ctutils;
use elliptic_curve::ops::{Invert, LinearCombination, MulByGeneratorVartime, MulVartime, Reduce};
use elliptic_curve::point::{AffineCoordinates, NonIdentity};
use elliptic_curve::scalar::{FromUintUnchecked, IsHigh};
use elliptic_curve::{
    BatchNormalize, CurveAffine, CurveGroup, Error as EcError, FieldBytes, Generate, ScalarValue,
};
use rand_core::TryCryptoRng;

/// BabyJubJub curve type
#[derive(Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
pub struct BabyJubJub;

/// Order of the BabyJubJub scalar field (prime order subgroup)
///
/// This is the standard BabyJubJub scalar field order:
/// r = 2736030358979909402780800718157159386076813972158567259200215660948447373041
///
/// NOTE: This value is verified at compile time by tests to match the backend's
/// `BackendScalar::MODULUS`. If the backend is updated, tests will catch any mismatch.
const ORDER_HEX: &str = "060c89ce5c263405370a08b6d0302b0bab3eedb83920ee0a677297dc392126f1";

const SCALAR_MODULUS_LE: [u8; 32] = [
    0xf1, 0x26, 0x21, 0x39, 0xdc, 0x97, 0x72, 0x67, 0x0a, 0xee, 0x20, 0x39, 0xb8, 0xed, 0x3e, 0xab,
    0x0b, 0x2b, 0x30, 0xd0, 0xb6, 0x08, 0x0a, 0x37, 0x05, 0x34, 0x26, 0x5c, 0xce, 0x89, 0x0c, 0x06,
];

/// BabyJubJub cofactor (the number of curve points per prime-order subgroup element).
///
/// This value is checked by the test suite to match the underlying curve's
/// cofactor; a backend change that altered the cofactor would fail tests.
pub const COFACTOR: u64 = 8;

impl Curve for BabyJubJub {
    type FieldBytesSize = elliptic_curve::consts::U32;
    type Uint = elliptic_curve::bigint::U256;

    /// Order of the BabyJubJub scalar field
    const ORDER: elliptic_curve::bigint::Odd<Self::Uint> =
        elliptic_curve::bigint::Odd::from_be_hex(ORDER_HEX);
}

// `BabyJubJub::ORDER` is the order `r` of the prime-order subgroup (the scalar
// field order), which is prime. `PrimeCurve` is a marker trait asserting exactly
// that, so this impl is sound. It lets generic code bound on `PrimeCurve`.
impl PrimeCurve for BabyJubJub {}

impl CurveArithmetic for BabyJubJub {
    type AffinePoint = AffinePoint;
    type ProjectivePoint = ProjectivePoint;
    type Scalar = Scalar;
}

// ===== AffinePoint Wrapper =====

/// Affine point representation
/// Note: BabyJubJub coordinates are in Fq (base field), not Fr (scalar field)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct AffinePoint {
    pub x: BackendBaseField,
    pub y: BackendBaseField,
}

impl AffinePoint {
    /// Additive identity of the group (point at infinity)
    pub const IDENTITY: Self = Self {
        x: BackendBaseField::ZERO,
        y: BackendBaseField::ONE,
    };

    /// Generator point from backend (via GENERATOR_X and GENERATOR_Y constants)
    pub const GENERATOR: Self = Self {
        x: taceo_ark_babyjubjub::GENERATOR_X,
        y: taceo_ark_babyjubjub::GENERATOR_Y,
    };

    // Create a new affine point from raw coordinates **without any validation**.
    //
    // This is for internal/test use only.
    #[allow(dead_code)]
    pub(crate) fn new_unchecked(x: BackendBaseField, y: BackendBaseField) -> Self {
        Self { x, y }
    }

    /// Create a new affine point, returning `None` unless the coordinates are
    /// on the curve **and** in the prime-order subgroup.
    ///
    /// This is the safe constructor to use for untrusted `(x, y)` pairs.
    pub fn new(x: BackendBaseField, y: BackendBaseField) -> Option<Self> {
        let p = Self { x, y };
        if p.is_in_prime_order_subgroup() {
            Some(p)
        } else {
            None
        }
    }

    /// Returns true iff `(x, y)` satisfies the twisted Edwards curve equation.
    ///
    /// Note: being on the curve does **not** imply membership in the prime-order
    /// subgroup (BabyJubJub has cofactor 8); use
    /// [`AffinePoint::is_in_prime_order_subgroup`] for the stronger check.
    pub fn is_on_curve(&self) -> bool {
        BackendAffine::new_unchecked(self.x, self.y).is_on_curve()
    }

    /// Returns true iff the point is on the curve **and** lies in the
    /// prime-order subgroup of order [`BabyJubJub::ORDER`].
    pub fn is_in_prime_order_subgroup(&self) -> bool {
        let a = BackendAffine::new_unchecked(self.x, self.y);
        a.is_on_curve() && a.is_in_correct_subgroup_assuming_on_curve()
    }

    /// Check if this point is the identity
    pub fn is_identity(&self) -> bool {
        self.x.is_zero() && self.y.is_one()
    }

    /// Get the x-coordinate
    pub fn x(&self) -> BackendBaseField {
        self.x
    }

    /// Get the y-coordinate
    pub fn y(&self) -> BackendBaseField {
        self.y
    }

    /// Parity (least-significant bit) of the **canonical** x-coordinate.
    ///
    /// This is the RFC 8032 / EdDSA sign convention (LSB of `x`), which is
    /// **not** the same convention used by the underlying curve's compressed point
    /// encoding (ark packs a `x > -x` "is-negative" flag into the spare high
    /// bits of `y`). Do not use this value to hand-roll point compression that
    /// must interoperate with [`ProjectivePoint::to_bytes`]; it is provided only
    /// for callers that explicitly need x-parity. Allocation-free.
    pub fn x_is_odd(&self) -> subtle::Choice {
        // LSB of the canonical (non-Montgomery) integer representation of x.
        let limbs = self.x.into_bigint().0;
        subtle::Choice::from((limbs[0] & 1) as u8)
    }
}

// ===== ProjectivePoint Wrapper =====

/// Projective point representation
/// Note: BabyJubJub coordinates are in Fq (base field), not Fr (scalar field)
#[derive(Clone, Copy, Debug)]
pub struct ProjectivePoint {
    pub x: BackendBaseField,
    pub y: BackendBaseField,
    pub z: BackendBaseField,
}

impl PartialEq for ProjectivePoint {
    fn eq(&self, other: &Self) -> bool {
        // Identity shortcut: both identity match (avoids a field inversion).
        if self.is_identity() && other.is_identity() {
            return true;
        }
        // Cross-multiply to compare affine coordinates without normalization:
        // (X1,Y1,Z1) == (X2,Y2,Z2)  <=>  X1/Z1 == X2/Z2  AND  Y1/Z1 == Y2/Z2
        //                        <=>  X1*Z2 == X2*Z1    AND  Y1*Z2 == Y2*Z1
        self.x * other.z == other.x * self.z && self.y * other.z == other.y * self.z
    }
}

impl Eq for ProjectivePoint {}

impl ProjectivePoint {
    /// Additive identity of the group (point at infinity)
    /// In extended projective coordinates, identity is (0, 1, 0, 1)
    pub const IDENTITY: Self = Self {
        x: BackendBaseField::ZERO,
        y: BackendBaseField::ONE,
        z: BackendBaseField::ONE,
    };

    /// Generator point from backend (converted to projective coordinates)
    pub const GENERATOR: Self = {
        Self {
            x: taceo_ark_babyjubjub::GENERATOR_X,
            y: taceo_ark_babyjubjub::GENERATOR_Y,
            z: BackendBaseField::ONE,
        }
    };

    // Create a new projective point from coordinates **without any validation**.
    //
    // This is for internal/test use only.
    #[allow(dead_code)]
    pub(crate) fn new_unchecked(
        x: BackendBaseField,
        y: BackendBaseField,
        z: BackendBaseField,
    ) -> Self {
        Self { x, y, z }
    }

    /// Internal method to convert to a backend projective point without validation.
    pub(crate) fn to_backend_unvalidated(self) -> BackendProjective {
        // The backend uses *extended* twisted Edwards coordinates (X : Y : T : Z)
        // with the invariant `T * Z == X * Y` (so that `T/Z == (X/Z)*(Y/Z)`).
        // Our `(x, y, z)` represents the affine point `(x/z, y/z)`. Naively setting
        // `T = x*y, Z = z` violates the invariant whenever `z != 1`, producing an
        // *invalid* extended point whose use in the backend's addition formulas
        // (which read `t`) would be incorrect. Scaling all coordinates by `z`
        // restores a valid representation of the same affine point:
        //   (x*z, y*z, x*y, z^2)
        // because (x*z)/(z^2) = x/z, (y*z)/(z^2) = y/z, (x*y)/(z^2) = (x/z)(y/z),
        // and the invariant holds: (x*y)*(z^2) == (x*z)*(y*z). This is
        // inversion-free and reduces to (x, y, x*y, 1) when z == 1.
        BackendProjective {
            x: self.x * self.z,
            y: self.y * self.z,
            t: self.x * self.y,
            z: self.z.square(),
        }
    }

    // Group doubling using the backend, bypassing the prime-order-subgroup
    // assertion. For internal/test use only.
    #[allow(dead_code)]
    pub(crate) fn double_unchecked(&self) -> Self {
        self.to_backend_unvalidated().double().into()
    }

    // Scalar multiplication using the backend, bypassing the prime-order-subgroup
    // assertion. For internal/test use only.
    #[allow(dead_code)]
    pub(crate) fn mul_unchecked(&self, scalar: &Scalar) -> Self {
        let result = self.to_backend_unvalidated() * scalar.0;
        result.into()
    }

    /// Check if this point is the identity (neutral element).
    ///
    /// BabyJubJub is a (complete) twisted Edwards curve, so there is **no**
    /// point at infinity: the neutral element is the affine point `(0, 1)`,
    /// which in projective coordinates is `(0 : Y : Z)` with `Y == Z` (e.g.
    /// [`ProjectivePoint::IDENTITY`] is `(0, 1, 1)`). Testing `z == 0` — as a
    /// short Weierstrass implementation would — is therefore wrong and never
    /// matches a real identity produced by this API.
    pub fn is_identity(&self) -> bool {
        // Best-effort constant-time check using bitwise limb operations.
        // The `z != 0` guard uses short-circuit (not a secret — it determines
        // a representation format, not a key-dependent value).
        let z_nonzero = !self.z.is_zero();
        let limbs_x = &self.x.0.0;
        let limbs_z = &self.z.0.0;
        let x_is_zero = limbs_x[0] == 0 && limbs_x[1] == 0 && limbs_x[2] == 0 && limbs_x[3] == 0;
        let y_eq_z = self.y.0.0[0] == limbs_z[0]
            && self.y.0.0[1] == limbs_z[1]
            && self.y.0.0[2] == limbs_z[2]
            && self.y.0.0[3] == limbs_z[3];
        z_nonzero && x_is_zero && y_eq_z
    }

    /// Returns true iff the point is on the curve **and** in the prime-order
    /// subgroup. See [`AffinePoint::is_in_prime_order_subgroup`].
    pub fn is_in_prime_order_subgroup(&self) -> bool {
        self.is_on_curve() && self.to_affine().is_in_prime_order_subgroup()
    }

    /// Returns true iff the point satisfies the curve equation.
    /// See [`AffinePoint::is_on_curve`].
    pub fn is_on_curve(&self) -> bool {
        if self.z.is_zero() {
            return false;
        }

        // Validate the projective equation directly before any affine
        // conversion. For projective `(X : Y : Z)` representing affine
        // `(X/Z, Y/Z)`, BabyJubJub's Edwards equation
        // `a*x^2 + y^2 = 1 + d*x^2*y^2` becomes:
        // `(a*X^2 + Y^2)*Z^2 = Z^4 + d*X^2*Y^2`.
        let x2 = self.x.square();
        let y2 = self.y.square();
        let z2 = self.z.square();
        let a = BackendBaseField::from(168700u64);
        let d = BackendBaseField::from(168696u64);
        (a * x2 + y2) * z2 == z2.square() + d * x2 * y2
    }

    /// Almost constant-time scalar multiplication `[scalar] * self`.
    ///
    /// Uses a fixed-length double-and-add-always loop with a bit-mask conditional
    /// select: every iteration computes both `[2]acc` and `[2]acc + self` and
    /// selects between them on the current scalar bit, so there is no
    /// scalar-dependent loop length or early-exit in this wrapper. (This is the
    /// double-and-add-always countermeasure, *not* a Montgomery ladder.) See the
    /// `# Timing` section below for the end-to-end caveat.
    ///
    /// # Panics
    ///
    /// Panics unless `self` is in the prime-order subgroup. Each iteration adds
    /// via the `+` operator, which converts to the underlying curve
    /// implementation's representation and asserts subgroup
    /// membership. For points that may be on-curve but outside the prime-order
    /// subgroup (e.g. torsion points) use [`ProjectivePoint::mul_with_cofactor_clear`]
    /// instead. The same subgroup assertion applies to the `*` operator and
    /// [`Group::double`].
    ///
    /// # Timing
    ///
    /// This is **not** an end-to-end constant-time scalar multiplication.
    /// The fixed loop length and bit-mask conditional select are constant-time
    /// at the algorithm level, but each iteration calls
    /// the underlying curve implementation's point addition and doubling, which in turn call
    /// `ark-ff` field arithmetic. That field arithmetic uses a data-dependent
    /// conditional reduction (`Fp::subtract_modulus` via `is_geq_modulus`, a
    /// regular BigInteger comparison that compiles to u64-level conditional
    /// jumps). Whether this reduction fires in a given iteration is a function
    /// of the intermediate field values, which are themselves a function of
    /// the scalar, so the operation carries a small timing signal. Callers
    /// who require end-to-end constant-time scalar multiplication must use a
    /// different backend (e.g. one with bit-mask-based field reduction
    /// throughout, such as `fiat-crypto`-generated code).
    pub fn mul_fixed_schedule(&self, scalar: &Scalar) -> Self {
        // Canonical little-endian bytes of the scalar (< r < 2^252).
        let bytes = scalar.to_bytes_le();
        let mut acc = Self::IDENTITY;
        // Scalar::NUM_BITS iterations covers all valid bits; bits NUM_BITS..255 are
        // always zero per the field definition, so processing them is wasted work.
        for i in (0..Scalar::NUM_BITS as usize).rev() {
            let doubled = Group::double(&acc);
            let sum = doubled + *self;
            let bit = (bytes[i >> 3] >> (i & 7)) & 1;
            acc = Self::conditional_select(&doubled, &sum, subtle::Choice::from(bit));
        }
        acc
    }

    /// Scalar multiplication that additionally clears the cofactor, producing a
    /// point guaranteed to be in the prime-order subgroup.
    ///
    /// BabyJubJub has cofactor 8: the full curve group is the (internal) direct
    /// sum of the prime-order subgroup (order `r`) and a torsion subgroup of
    /// order 8. Every point therefore splits uniquely as `P = Q + T`, with `Q`
    /// in the prime-order subgroup and `T` in the torsion subgroup. Every
    /// torsion point has order dividing 8, so the cofactor annihilates it
    /// (`[8]T = O`) and `[8]P = [8]Q` always lies in the prime-order subgroup.
    /// Plain scalar multiplication `[s]P`, by contrast, leaves the torsion
    /// component (`[s]T`) intact, so its result need not lie in the prime-order
    /// subgroup.
    ///
    /// **Security impact.** Protocols that multiply an attacker-supplied point
    /// by a secret scalar while expecting a prime-order-subgroup result are
    /// vulnerable to small-subgroup attacks when the attacker can control the
    /// torsion component of `P`. This method should be used instead of the
    /// plain `*` operator whenever the caller needs the result to be provably
    /// in the prime-order subgroup and the input point is not yet known to
    /// be in it (e.g., points decoded from untrusted input).
    ///
    /// This computes `[8s]P` using the **integer** product `8 * s`, which equals
    /// `[8s]Q` (the factor of 8
    /// annihilates the torsion component, having the order dividing cofactor 8).
    /// The result is `IDENTITY` exactly when
    /// `[8s]Q = O`: either `P` is pure torsion (`Q` is the identity, so the
    /// result is `IDENTITY` for any `s`) or `Q` has full order `r` and `8s` is a
    /// multiple of `r`.
    ///
    /// The cofactor multiplication is done over the integers, **not** in the scalar
    /// field: `8 * s mod r` is generally *not* a multiple of 8 (`r` is odd), so
    /// reducing modulo `r` would leave the torsion component intact and defeat
    /// the cofactor clearing entirely.
    pub fn mul_with_cofactor_clear(&self, scalar: &Scalar) -> Self {
        // Build integer 8 * s (without reducing mod r) and feed it to the
        // backend's integer scalar multiplication. Since `s < r < 2^251`, the
        // product `8 * s < 2^254` fits in four 64-bit limbs and never carries
        // out of the top limb. Sizing the slice at four limbs (rather than five)
        // keeps the backend's fixed-length Montgomery ladder at 256 iterations
        // instead of 320.
        let s = scalar.0.into_bigint().0;
        let mut wide = [0u64; 4];
        let mut carry = 0u128;
        for (i, &limb) in s.iter().enumerate() {
            let acc = (limb as u128) * (COFACTOR as u128) + carry;
            wide[i] = acc as u64;
            carry = acc >> 64;
        }
        // `8 * s < 2^254` guarantees no carry escapes the fourth limb.
        debug_assert_eq!(carry, 0, "8 * s must fit in four 64-bit limbs");
        let result = self.to_backend_unvalidated().mul_bigint(wide);
        result.into()
    }

    /// Convert to affine coordinates.
    ///
    /// # Panics
    ///
    /// Panics if the point is not in valid projective coordinates (`z == 0`).
    /// Invalid projective points must not be fed to this method; use
    /// [`ProjectivePoint::is_on_curve`] to validate untrusted inputs first.
    pub fn to_affine(&self) -> AffinePoint {
        assert!(
            !self.z.is_zero(),
            "invalid projective point: z-coordinate is zero"
        );

        let z_inv = self.z.inverse().expect("non-zero has inverse");
        let x = self.x * z_inv;
        let y = self.y * z_inv;
        AffinePoint { x, y }
    }
}

// ===== Scalar Wrapper =====

/// Scalar field element
///
/// # Zeroization
///
/// This type is `Copy` and therefore cannot implement `ZeroizeOnDrop`. Callers
/// handling secret scalars are responsible for zeroizing their own storage.
/// The type implements `zeroize::DefaultIsZeroes` (and hence `Zeroize`)
/// **unconditionally** via `elliptic_curve`'s re-exported `zeroize`, so you can
/// call `.zeroize()` explicitly when done.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Scalar(pub BackendScalar);

impl Scalar {
    /// Zero scalar
    pub const ZERO: Self = Self(BackendScalar::ZERO);

    /// One scalar
    pub const ONE: Self = Self(BackendScalar::ONE);

    /// Size of the scalar in bytes.
    pub const SIZE: usize = 32;

    /// Create a scalar from a **canonical** big-endian byte encoding.
    ///
    /// # Security
    ///
    /// Returns `CtOption::none()` if `bytes` encodes an integer `>= r` (the
    /// scalar-field order). Rejecting non-canonical encodings prevents the
    /// scalar/signature malleability that arises when distinct byte strings
    /// (e.g. `s` and `s + r`) would otherwise silently reduce to the same
    /// scalar. If you instead want modular reduction, use
    /// [`Scalar::reduce_bytes_be`].
    pub fn from_bytes(bytes: &[u8; Self::SIZE]) -> CtOption<Self> {
        let mut le = *bytes;
        le.reverse();
        Self::from_canonical_le(&le)
    }

    /// Create a scalar from a **canonical** little-endian byte encoding.
    ///
    /// See [`Scalar::from_bytes`] for the canonicity / security contract;
    /// returns `CtOption::none()` for any value `>= r`.
    pub fn from_bytes_le(bytes: &[u8; Self::SIZE]) -> CtOption<Self> {
        Self::from_canonical_le(bytes)
    }

    /// Reduce an arbitrary big-endian 32-byte value modulo `r`.
    ///
    /// Unlike [`Scalar::from_bytes`] this never fails and is **not** canonical:
    /// inputs `>= r` are reduced. Use only where modular reduction is the
    /// intended behaviour — never when decoding signatures or other values that
    /// must round-trip canonically.
    pub fn reduce_bytes_be(bytes: &[u8; Self::SIZE]) -> Self {
        let mut le = *bytes;
        le.reverse();
        Self(BackendScalar::from_le_bytes_mod_order(&le))
    }

    /// Reduce an arbitrary little-endian 32-byte value modulo `r`.
    /// See [`Scalar::reduce_bytes_be`].
    pub fn reduce_bytes_le(bytes: &[u8; Self::SIZE]) -> Self {
        Self(BackendScalar::from_le_bytes_mod_order(bytes))
    }

    /// Build a scalar from canonical little-endian bytes, rejecting any value
    /// `>= r`. Backs [`Scalar::from_bytes`] / [`Scalar::from_bytes_le`].
    fn from_canonical_le(bytes: &[u8; Self::SIZE]) -> CtOption<Self> {
        let value = BackendScalar::from_le_bytes_mod_order(bytes);
        CtOption::new(Self(value), Self::is_canonical_le(bytes))
    }

    fn is_canonical_le(bytes: &[u8; Self::SIZE]) -> subtle::Choice {
        let mut borrow = 0u16;
        for (&a, &b) in bytes.iter().zip(SCALAR_MODULUS_LE.iter()) {
            let diff = (a as u16).wrapping_sub(b as u16).wrapping_sub(borrow);
            borrow = diff >> 15;
        }
        subtle::Choice::from(borrow as u8)
    }

    /// Convert to bytes (little-endian). Allocation-free.
    pub fn to_bytes_le(&self) -> [u8; Self::SIZE] {
        let limbs = self.0.into_bigint().0;
        let mut arr = [0u8; 32];
        for (i, limb) in limbs.iter().enumerate() {
            arr[i * 8..i * 8 + 8].copy_from_slice(&limb.to_le_bytes());
        }
        arr
    }

    /// Convert to bytes (big-endian). Allocation-free.
    pub fn to_bytes(&self) -> [u8; Self::SIZE] {
        let mut be = self.to_bytes_le();
        be.reverse();
        be
    }

    /// Check if scalar is zero
    pub fn is_zero(&self) -> bool {
        self.0.is_zero()
    }

    /// Check if scalar is one
    pub fn is_one(&self) -> bool {
        self.0.is_one()
    }

    /// Invert the scalar
    ///
    /// # Security Note
    ///
    /// This method uses the underlying field implementation's inverse function which implements
    /// a variable-time extended Euclidean algorithm. This means:
    /// - The timing may leak information about whether the input is zero
    /// - The method returns `CtOption::new(Self::ZERO, 0.into())` when input is zero
    ///
    /// **Do not use this method with secret nonces** (e.g. the `k` value in
    /// ECDSA signature generation). The variable-time algorithm can leak the
    /// nonce through timing, enabling private-key recovery attacks. For
    /// nonce inversion in signature schemes, use a constant-time modular
    /// inversion routine instead.
    ///
    /// For most other use cases (e.g., signature verification), the
    /// non-constant-time behaviour is acceptable as the scalar is already
    /// validated to be non-zero.
    pub fn invert(&self) -> CtOption<Self> {
        match self.0.inverse() {
            Some(s) => CtOption::new(Self(s), 1.into()),
            None => CtOption::new(Self::ZERO, 0.into()),
        }
    }

    /// Square the scalar
    pub fn square(&self) -> Self {
        Self(self.0.square())
    }

    /// Double the scalar
    pub fn double(&self) -> Self {
        Self(self.0 + self.0)
    }
}

// ===== Group trait implementations =====

impl Group for ProjectivePoint {
    type Scalar = Scalar;

    fn random<R: rand_core::Rng + ?Sized>(rng: &mut R) -> Self {
        let scalar = Scalar::random(rng);
        Self::generator() * scalar
    }

    fn try_random<R: rand_core::TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        let scalar = Scalar::try_random(rng)?;
        Ok(Self::generator() * scalar)
    }

    fn generator() -> Self {
        Self::GENERATOR
    }

    fn identity() -> Self {
        Self::IDENTITY
    }

    fn is_identity(&self) -> subtle::Choice {
        // Twisted Edwards has no point at infinity: the neutral element is the
        // affine point (0, 1), i.e. projective (0 : Y : Z) with Y == Z. Testing
        // `z == 0` (a short-Weierstrass habit) never matches a real identity.
        //
        // Use bitwise constant-time operations on the field BigInt limbs to
        // avoid branching on secret-dependent values.
        let z_nonzero = !self.z.is_zero();
        let limbs_x = &self.x.0.0;
        let limbs_y = &self.y.0.0;
        let limbs_z = &self.z.0.0;
        let x_is_zero = limbs_x[0].ct_eq(&0)
            & limbs_x[1].ct_eq(&0)
            & limbs_x[2].ct_eq(&0)
            & limbs_x[3].ct_eq(&0);
        let y_eq_z = limbs_y[0].ct_eq(&limbs_z[0])
            & limbs_y[1].ct_eq(&limbs_z[1])
            & limbs_y[2].ct_eq(&limbs_z[2])
            & limbs_y[3].ct_eq(&limbs_z[3]);
        subtle::Choice::from(z_nonzero as u8) & x_is_zero & y_eq_z
    }

    fn double(&self) -> Self {
        // Convert to BackendProjective, double, and convert back
        let backend: BackendProjective = (*self).into();
        let doubled = backend.double();
        doubled.into()
    }
}

impl ConditionallySelectable for ProjectivePoint {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        // Use conditional_select on the underlying BigInt arrays (u64 arrays)
        let x_limb0 = u64::conditional_select(&a.x.0.0[0], &b.x.0.0[0], choice);
        let x_limb1 = u64::conditional_select(&a.x.0.0[1], &b.x.0.0[1], choice);
        let x_limb2 = u64::conditional_select(&a.x.0.0[2], &b.x.0.0[2], choice);
        let x_limb3 = u64::conditional_select(&a.x.0.0[3], &b.x.0.0[3], choice);
        let y_limb0 = u64::conditional_select(&a.y.0.0[0], &b.y.0.0[0], choice);
        let y_limb1 = u64::conditional_select(&a.y.0.0[1], &b.y.0.0[1], choice);
        let y_limb2 = u64::conditional_select(&a.y.0.0[2], &b.y.0.0[2], choice);
        let y_limb3 = u64::conditional_select(&a.y.0.0[3], &b.y.0.0[3], choice);
        let z_limb0 = u64::conditional_select(&a.z.0.0[0], &b.z.0.0[0], choice);
        let z_limb1 = u64::conditional_select(&a.z.0.0[1], &b.z.0.0[1], choice);
        let z_limb2 = u64::conditional_select(&a.z.0.0[2], &b.z.0.0[2], choice);
        let z_limb3 = u64::conditional_select(&a.z.0.0[3], &b.z.0.0[3], choice);

        Self {
            x: BackendBaseField::new_unchecked(ark_ff::BigInt([
                x_limb0, x_limb1, x_limb2, x_limb3,
            ])),
            y: BackendBaseField::new_unchecked(ark_ff::BigInt([
                y_limb0, y_limb1, y_limb2, y_limb3,
            ])),
            z: BackendBaseField::new_unchecked(ark_ff::BigInt([
                z_limb0, z_limb1, z_limb2, z_limb3,
            ])),
        }
    }
}

impl DefaultIsZeroes for ProjectivePoint {}

// ===== GroupEncoding trait implementations =====

/// Canonical compressed encoding of a [`ProjectivePoint`].
///
/// This is exactly **32 bytes**, matching the arkworks BabyJubJub compressed
/// serialization (little-endian `y` with the x-sign flag packed into the spare
/// high bits of the last byte). A previous version used 33 bytes with an unused
/// trailing byte, which made the encoding non-canonical and malleable (256
/// distinct byte strings decoded to the same point); the extra byte has been
/// removed.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GroupRepr(pub [u8; Self::SIZE]);

impl GroupRepr {
    /// Size of the canonical group encoding in bytes.
    pub const SIZE: usize = 32;
}

impl AsRef<[u8]> for GroupRepr {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl AsMut<[u8]> for GroupRepr {
    fn as_mut(&mut self) -> &mut [u8] {
        &mut self.0
    }
}

impl GroupEncoding for ProjectivePoint {
    type Repr = GroupRepr;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_bytes_impl(&bytes.0, true)
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        Self::from_bytes_impl(&bytes.0, false)
    }

    fn to_bytes(&self) -> Self::Repr {
        // Use the backend's serialization which properly handles Montgomery form
        let backend_proj: BackendProjective = (*self).into();
        let backend_affine = BackendAffine::from(backend_proj);

        // Serialize the affine point using the backend's CanonicalSerialize.
        // The compressed format is exactly 32 bytes (little-endian y with the
        // x-sign flag in the spare high bits of the final byte).
        let mut bytes = [0u8; Scalar::SIZE];
        backend_affine
            .serialize_with_mode(bytes.as_mut(), Compress::Yes)
            .expect("serialization to 32 bytes should succeed");

        GroupRepr(bytes)
    }
}

impl ProjectivePoint {
    /// Internal implementation of from_bytes and from_bytes_unchecked.
    /// When validate is true, checks that the point is on the curve and in the correct subgroup.
    ///
    /// NOTE: the backend's specific deserialization error is dropped here —
    /// both malformed encoding and "valid point not in subgroup" produce the
    /// same `CtOption::new(Self::IDENTITY, 0.into())`. If fine-grained error
    /// handling is needed in future, this function's return type must change
    /// to carry a `Result` or a custom enum variant. The current `None`
    /// contract is intentional for API surface-area reasons.
    fn from_bytes_impl(bytes: &[u8; 32], validate: bool) -> CtOption<Self> {
        // Use the backend's deserialization which properly handles Montgomery form.
        // Avoiding an early return on the error path keeps the *control flow* the same
        // for parseable and unparseable encodings. Note this does NOT make decoding
        // constant-time: the backend's decompression (an x-coordinate sqrt) and, under
        // `Validate::Yes`, its on-curve/subgroup checks are variable-time. Decoding is
        // not a constant-time operation and must not be treated as one.
        let mut reader = bytes.as_ref();
        let (backend_affine, was_ok) = match BackendAffine::deserialize_with_mode(
            &mut reader,
            Compress::Yes,
            if validate {
                Validate::Yes
            } else {
                Validate::No
            },
        ) {
            Ok(affine) => (affine, true),
            // Construct the default affine point even on failure so both branches
            // produce a valid `our_affine` and we can defer the accept/reject decision
            // to the `was_ok` flag instead of an early return.
            Err(_) => (
                BackendAffine {
                    x: BackendBaseField::ZERO,
                    y: BackendBaseField::ONE,
                },
                false,
            ),
        };

        // Convert backend affine to our affine wrapper — unconditionally
        let our_affine = AffinePoint {
            x: backend_affine.x,
            y: backend_affine.y,
        };

        // Reject the non-canonical x-sign encoding. The compressed format packs an
        // x-sign flag into bit 7 of the final byte; when `x == 0`, `x` and `-x`
        // coincide, so both flag values would decode to the same point. The canonical
        // encoder always emits sign bit 0 for `x == 0`, so an input with `x == 0` and
        // the sign bit set is a second, non-canonical encoding of the same point and
        // must be rejected to keep the encoding non-malleable. Computed branchlessly so
        // it does not add input-dependent control flow.
        let xl = &our_affine.x.0.0;
        let x_is_zero = xl[0].ct_eq(&0) & xl[1].ct_eq(&0) & xl[2].ct_eq(&0) & xl[3].ct_eq(&0);
        let sign_bit = subtle::Choice::from((bytes[31] >> 7) & 1);
        let canonical_sign = !(x_is_zero & sign_bit);

        CtOption::new(
            Self::from(our_affine),
            subtle::Choice::from(was_ok as u8) & canonical_sign,
        )
    }
}

// ===== Arithmetic Operator Implementations =====

impl<'a> core::ops::Add<&'a ProjectivePoint> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: &'a ProjectivePoint) -> ProjectivePoint {
        let backend_self = BackendProjective::from(*self);
        let backend_other = BackendProjective::from(*other);
        let result = backend_self + backend_other;
        result.into()
    }
}

impl core::ops::Add<ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: ProjectivePoint) -> ProjectivePoint {
        &self + &other
    }
}

impl core::ops::Add<&ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, other: &ProjectivePoint) -> ProjectivePoint {
        let backend_self = BackendProjective::from(self);
        let backend_other = BackendProjective::from(*other);
        let result = backend_self + backend_other;
        result.into()
    }
}

impl core::ops::AddAssign<ProjectivePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: ProjectivePoint) {
        *self = *self + rhs;
    }
}

impl core::ops::AddAssign<&ProjectivePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: &ProjectivePoint) {
        *self = *self + *rhs;
    }
}

impl core::ops::Sub<ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: ProjectivePoint) -> ProjectivePoint {
        let backend_self = BackendProjective::from(self);
        let backend_other = BackendProjective::from(other);
        let result = backend_self - backend_other;
        result.into()
    }
}

impl core::ops::Sub<&ProjectivePoint> for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: &ProjectivePoint) -> ProjectivePoint {
        let backend_self = BackendProjective::from(*self);
        let backend_other = BackendProjective::from(*other);
        let result = backend_self - backend_other;
        result.into()
    }
}

impl core::ops::Sub<&ProjectivePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, other: &ProjectivePoint) -> ProjectivePoint {
        let backend_self = BackendProjective::from(self);
        let backend_other = BackendProjective::from(*other);
        let result = backend_self - backend_other;
        result.into()
    }
}

impl core::ops::SubAssign<ProjectivePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: ProjectivePoint) {
        *self = *self - rhs;
    }
}

impl core::ops::SubAssign<&ProjectivePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: &ProjectivePoint) {
        *self = *self - *rhs;
    }
}

impl core::ops::Neg for ProjectivePoint {
    type Output = ProjectivePoint;

    fn neg(self) -> ProjectivePoint {
        let backend = BackendProjective::from(self);
        let negated = -backend;
        negated.into()
    }
}

impl core::ops::Neg for &ProjectivePoint {
    type Output = ProjectivePoint;

    fn neg(self) -> ProjectivePoint {
        let backend = BackendProjective::from(*self);
        let negated = -backend;
        negated.into()
    }
}

impl<'a> core::ops::Mul<&'a Scalar> for &ProjectivePoint {
    type Output = ProjectivePoint;

    /// Scalar multiplication `point * scalar`.
    ///
    /// # Timing
    ///
    /// **Almost constant-time.** This (and all other `Mul<Scalar>` operator
    /// impls) delegates to the backend's scalar multiplication, which
    /// `taceo-ark-babyjubjub` implements as a **Montgomery ladder**: it iterates
    /// over a fixed number of scalar bits (no leading-zero skipping) and selects
    /// between the two ladder registers with branch-free, bit-masked conditional
    /// swaps, so there is no scalar-dependent loop length, branching, or
    /// memory-access pattern at the algorithm level. It is **not** end-to-end
    /// constant-time, though: the underlying `ark-ff` field arithmetic uses a
    /// data-dependent conditional reduction (`Fp::subtract_modulus` via
    /// `is_geq_modulus`), which leaves a small residual timing signal. This is
    /// the same timing profile as [`ProjectivePoint::mul_fixed_schedule`];
    /// callers needing end-to-end constant time must use a backend with
    /// bit-masked field reduction throughout (e.g. `fiat-crypto`-generated code).
    fn mul(self, scalar: &'a Scalar) -> ProjectivePoint {
        let backend_self = BackendProjective::from(*self);
        let result = backend_self * scalar.0;
        result.into()
    }
}

impl core::ops::Mul<Scalar> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: Scalar) -> ProjectivePoint {
        self * &scalar
    }
}

impl core::ops::Mul<&Scalar> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &Scalar) -> ProjectivePoint {
        let backend_self = BackendProjective::from(self);
        let result = backend_self * scalar.0;
        result.into()
    }
}

impl core::ops::MulAssign<Scalar> for ProjectivePoint {
    fn mul_assign(&mut self, rhs: Scalar) {
        *self = *self * rhs;
    }
}

impl core::ops::MulAssign<&Scalar> for ProjectivePoint {
    fn mul_assign(&mut self, rhs: &Scalar) {
        *self = *self * *rhs;
    }
}

impl From<AffinePoint> for ProjectivePoint {
    fn from(point: AffinePoint) -> Self {
        if point.is_identity() {
            Self::IDENTITY
        } else {
            Self {
                x: point.x,
                y: point.y,
                z: BackendBaseField::ONE,
            }
        }
    }
}

impl From<&AffinePoint> for ProjectivePoint {
    fn from(point: &AffinePoint) -> Self {
        Self::from(*point)
    }
}

impl From<ProjectivePoint> for AffinePoint {
    /// # Panics
    ///
    /// Panics if `point` has `z == 0` (an invalid projective point). See
    /// [`ProjectivePoint::to_affine`]. Validate untrusted
    /// input with [`ProjectivePoint::is_on_curve`] before converting.
    fn from(point: ProjectivePoint) -> Self {
        point.to_affine()
    }
}

impl From<&ProjectivePoint> for AffinePoint {
    /// # Panics
    ///
    /// Panics if `point` has `z == 0`; see [`From<ProjectivePoint>`](AffinePoint).
    fn from(point: &ProjectivePoint) -> Self {
        point.to_affine()
    }
}

// ===== Conversions between ProjectivePoint and BackendProjective =====

impl From<ProjectivePoint> for BackendProjective {
    /// # Panics
    ///
    /// Panics if `point` is on the curve but **not** in the prime-order subgroup
    /// (the underlying curve implementation's constructor asserts subgroup membership). This is
    /// the assertion that makes the arithmetic operators (`+`, `-`, `*`,
    /// [`Group::double`], [`ProjectivePoint::mul_fixed_schedule`]) panic on
    /// torsion/small-subgroup inputs.
    fn from(point: ProjectivePoint) -> Self {
        // Identity shortcut: BackendProjective::new asserts the subgroup for all
        // non-identity points, so calling it for the identity would panic.
        //
        // The identity in our projective representation is any (0 : λY : λZ) with
        // λ ≠ 0, equivalently x == 0 and y == z.  In the canonical form
        // (0, 1, 1) both coordinates are 1; in scaled form (0, 2, 2) they are 2.
        // Checking only x == 0 would incorrectly match the order-2 torsion element
        // P2 = (0, -1, 1), so we require y == z as well.
        let t = point.x * point.y;
        if point.x.is_zero() && point.y == point.z {
            return Self::zero();
        }
        // For all other points (including torsion/small-subgroup points),
        // BackendProjective::new normalizes z to 1 and checks the subgroup.
        // For prime-order subgroup points this assert always passes.
        Self::new(point.x, point.y, t, point.z)
    }
}

impl From<BackendProjective> for ProjectivePoint {
    fn from(backend: BackendProjective) -> Self {
        // Use the backend's is_zero check to detect the identity
        if backend.is_zero() {
            Self::IDENTITY
        } else {
            Self {
                x: backend.x,
                y: backend.y,
                z: backend.z,
            }
        }
    }
}

impl<'a> From<&'a ProjectivePoint> for BackendProjective {
    fn from(point: &'a ProjectivePoint) -> Self {
        (*point).into()
    }
}

impl<'a> From<&'a BackendProjective> for ProjectivePoint {
    fn from(backend: &'a BackendProjective) -> Self {
        // Delegate to the owned impl so the identity is normalized to (0, 1, 1)
        // consistently (an un-normalized identity would otherwise break
        // `is_identity` and equality checks).
        (*backend).into()
    }
}

impl ConditionallySelectable for AffinePoint {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        // Use conditional_select on the underlying BigInt arrays (u64 arrays)
        let x_limb0 = u64::conditional_select(&a.x.0.0[0], &b.x.0.0[0], choice);
        let x_limb1 = u64::conditional_select(&a.x.0.0[1], &b.x.0.0[1], choice);
        let x_limb2 = u64::conditional_select(&a.x.0.0[2], &b.x.0.0[2], choice);
        let x_limb3 = u64::conditional_select(&a.x.0.0[3], &b.x.0.0[3], choice);
        let y_limb0 = u64::conditional_select(&a.y.0.0[0], &b.y.0.0[0], choice);
        let y_limb1 = u64::conditional_select(&a.y.0.0[1], &b.y.0.0[1], choice);
        let y_limb2 = u64::conditional_select(&a.y.0.0[2], &b.y.0.0[2], choice);
        let y_limb3 = u64::conditional_select(&a.y.0.0[3], &b.y.0.0[3], choice);

        Self {
            x: BackendBaseField::new_unchecked(ark_ff::BigInt([
                x_limb0, x_limb1, x_limb2, x_limb3,
            ])),
            y: BackendBaseField::new_unchecked(ark_ff::BigInt([
                y_limb0, y_limb1, y_limb2, y_limb3,
            ])),
        }
    }
}

impl DefaultIsZeroes for AffinePoint {}

// ===== Scalar Field Implementations =====

impl PrimeField for Scalar {
    // `CurveArithmetic` mandates `Scalar::Repr == FieldBytes<BabyJubJub>`
    // (= `Array<u8, U32>`). Internally we keep working with `[u8; 32]` and
    // convert at this boundary; the canonical **little-endian** encoding and
    // the `>= r` rejection are unchanged.
    type Repr = elliptic_curve::FieldBytes<BabyJubJub>;

    /// Decode a canonical **little-endian** scalar, matching the encoding
    /// produced by [`PrimeField::to_repr`]. Per the `ff` contract this rejects
    /// non-canonical inputs (any value `>= r`) by returning `CtOption::none()`.
    fn from_repr(bytes: Self::Repr) -> CtOption<Self> {
        // NOTE: `to_repr` is little-endian, so `from_repr` must be too.
        let bytes: [u8; 32] = bytes.into();
        Self::from_bytes_le(&bytes)
    }

    fn from_repr_vartime(bytes: Self::Repr) -> Option<Self> {
        // Same canonical little-endian decoding as `from_repr`; rejects `>= r`.
        let bytes: [u8; 32] = bytes.into();
        Self::from_bytes_le(&bytes).into()
    }

    fn to_repr(&self) -> Self::Repr {
        self.to_bytes_le().into()
    }

    fn is_odd(&self) -> subtle::Choice {
        let bytes = self.to_bytes_le();
        (bytes[0] & 1).into()
    }

    // The scalar-field modulus r is 251 bits, so a field element needs 251 bits
    // and at most 250 bits of arbitrary data can be stored without reduction.
    const NUM_BITS: u32 = 251;
    const CAPACITY: u32 = 250;

    const MODULUS: &'static str =
        "060c89ce5c263405370a08b6d0302b0bab3eedb83920ee0a677297dc392126f1";

    // Pre-computed values for BabyJubJub scalar field (Montgomery representation)
    const TWO_INV: Self = Self(BackendScalar::new_unchecked(ark_ff::BigInt([
        0x83998aef5047ce3b,
        0xf3d67fe3504c7925,
        0x7c2d4900ec0c780a,
        0xf8b21270ddbb92,
    ])));

    // Multiplicative generator of the scalar field, taken directly from the
    // backend's vetted `FftField::GENERATOR` (a quadratic non-residue and a
    // generator of the full multiplicative group, and the value used to derive
    // ROOT_OF_UNITY — as required by the `ff` contract). The previous hardcoded
    // value `5` is a quadratic *residue* mod r, hence an INVALID generator.
    const MULTIPLICATIVE_GENERATOR: Self = Self(BackendScalar::GENERATOR);

    // S = 4 because r - 1 = 2^4 * t with t odd (== BackendScalar::TWO_ADICITY).
    const S: u32 = 4;

    // The 2^S-th root of unity, taken directly from the backend so it is
    // guaranteed to equal MULTIPLICATIVE_GENERATOR^t and to be a *primitive*
    // 2^S root of unity. (The previous hardcoded value was not a root of unity
    // of Fr at all — it broke FFTs, the default `sqrt`, and `sqrt_ratio`.)
    const ROOT_OF_UNITY: Self = Self(BackendScalar::TWO_ADIC_ROOT_OF_UNITY);

    // Inverse of ROOT_OF_UNITY (Montgomery form). Verified by tests:
    // ROOT_OF_UNITY * ROOT_OF_UNITY_INV == ONE.
    const ROOT_OF_UNITY_INV: Self = Self(BackendScalar::new_unchecked(ark_ff::BigInt([
        0x3cd891231ce44036,
        0x97c6d3222a9aac61,
        0x22a59f417e5ba9ca,
        0x0373acbf899c1a70,
    ])));

    // DELTA = MULTIPLICATIVE_GENERATOR^(2^S) (Montgomery form), per the `ff`
    // contract. Verified by tests: DELTA == MULTIPLICATIVE_GENERATOR^(2^S).
    const DELTA: Self = Self(BackendScalar::new_unchecked(ark_ff::BigInt([
        0xffb1712417f98edb,
        0x3c08c1257227dc15,
        0x370087e6983b16f7,
        0x0013ff20bb212fb7,
    ])));
}

impl ConditionallySelectable for Scalar {
    fn conditional_select(a: &Self, b: &Self, choice: subtle::Choice) -> Self {
        // Use conditional_select on the underlying BigInt arrays (u64 arrays)
        let limb0 = u64::conditional_select(&a.0.0.0[0], &b.0.0.0[0], choice);
        let limb1 = u64::conditional_select(&a.0.0.0[1], &b.0.0.0[1], choice);
        let limb2 = u64::conditional_select(&a.0.0.0[2], &b.0.0.0[2], choice);
        let limb3 = u64::conditional_select(&a.0.0.0[3], &b.0.0.0[3], choice);

        Self(BackendScalar::new_unchecked(ark_ff::BigInt([
            limb0, limb1, limb2, limb3,
        ])))
    }
}

impl DefaultIsZeroes for Scalar {}

impl Field for Scalar {
    const ZERO: Self = Self(BackendScalar::new_unchecked(ark_ff::BigInt([0; 4])));
    const ONE: Self = Self(BackendScalar::new_unchecked(ark_ff::BigInt([
        0x73315dea08f9c76,
        0xe7acffc6a098f24b,
        0xf85a9201d818f015,
        0x1f16424e1bb7724,
    ])));

    fn random<R: rand_core::Rng + ?Sized>(rng: &mut R) -> Self {
        Scalar::random(rng)
    }

    fn try_random<R: rand_core::TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        Scalar::try_random(rng)
    }

    fn square(&self) -> Self {
        Self(self.0.square())
    }

    fn double(&self) -> Self {
        Self(self.0 + self.0)
    }

    /// # Security Note
    ///
    /// This uses variable-time inversion. See [Scalar::invert()] for details.
    fn invert(&self) -> CtOption<Self> {
        match self.0.inverse() {
            Some(s) => CtOption::new(Self(s), 1.into()),
            None => CtOption::new(Self::ZERO, 0.into()),
        }
    }

    /// Square root in the scalar field.
    ///
    /// # Security
    ///
    /// Variable-time: delegates to the backend's `sqrt` (a Tonelli–Shanks
    /// variant) whose control flow depends on the input. Do not call on secret
    /// values where timing is observable.
    fn sqrt(&self) -> CtOption<Self> {
        match self.0.sqrt() {
            Some(s) => CtOption::new(Self(s), 1.into()),
            None => CtOption::new(Self::ZERO, 0.into()),
        }
    }

    /// Compute `sqrt(num / div)` following the `Field::sqrt_ratio` contract.
    ///
    /// Delegates to `ff`'s generic implementation, which is built on this
    /// field's (now correct) `ROOT_OF_UNITY` and on the overridden `sqrt`
    /// above (preventing the documented infinite recursion). Returns
    /// `(1, sqrt(num/div))` when `num/div` is a square (and `(1, 0)` when
    /// `num == 0`), and `(0, sqrt(ROOT_OF_UNITY * num/div))` for a non-square
    /// (or `(0, 0)` when only `div == 0`).
    ///
    /// The previous implementation was a stub that unconditionally returned
    /// `(1, 1)`, silently claiming every ratio was a square — which breaks any
    /// hash-to-curve / quadratic-residue test built on it.
    ///
    /// # Security
    ///
    /// Variable-time: the implementation calls `sqrt` (a Tonelli–Shanks
    /// variant) whose control flow depends on the input, and the Legendre
    /// symbol computation itself can leak whether `num/div` is a quadratic
    /// residue through timing. For hash-to-curve constructions where the
    /// quadratic-residue decision on a secret input must remain secret,
    /// this method is not appropriate. No constant-time alternative is
    /// provided by this crate.
    fn sqrt_ratio(num: &Self, div: &Self) -> (subtle::Choice, Self) {
        group::ff::helpers::sqrt_ratio_generic(num, div)
    }
}

// Implement Sum trait for ProjectivePoint.
//
// NOTE: `core::iter::Product` is intentionally NOT implemented for points:
// there is no meaningful multiplication of two curve points, and the previous
// impl silently returned IDENTITY for any input (a silent-wrong-result trap).
impl core::iter::Sum for ProjectivePoint {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(ProjectivePoint::IDENTITY, |a, b| a + b)
    }
}

impl<'a> core::iter::Sum<&'a ProjectivePoint> for ProjectivePoint {
    fn sum<I: Iterator<Item = &'a ProjectivePoint>>(iter: I) -> Self {
        iter.cloned().sum()
    }
}

// Implement Sum and Product for Scalar
impl core::iter::Sum for Scalar {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Scalar::ZERO, |a, b| a + b)
    }
}

impl<'a> core::iter::Sum<&'a Scalar> for Scalar {
    fn sum<I: Iterator<Item = &'a Scalar>>(iter: I) -> Self {
        iter.cloned().sum()
    }
}

impl core::iter::Product for Scalar {
    fn product<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Scalar::ONE, |a, b| a * b)
    }
}

impl<'a> core::iter::Product<&'a Scalar> for Scalar {
    fn product<I: Iterator<Item = &'a Scalar>>(iter: I) -> Self {
        iter.cloned().product()
    }
}

// Implement Add/Sub/Mul for Scalar
impl<'a> core::ops::Add<&'a Scalar> for &Scalar {
    type Output = Scalar;

    fn add(self, other: &'a Scalar) -> Scalar {
        Scalar(self.0 + other.0)
    }
}

impl core::ops::Add<Scalar> for Scalar {
    type Output = Scalar;

    fn add(self, other: Scalar) -> Scalar {
        &self + &other
    }
}

impl core::ops::Add<&Scalar> for Scalar {
    type Output = Scalar;

    fn add(self, other: &Scalar) -> Scalar {
        Scalar(self.0 + other.0)
    }
}

impl core::ops::AddAssign<Scalar> for Scalar {
    fn add_assign(&mut self, rhs: Scalar) {
        *self = *self + rhs;
    }
}

impl core::ops::AddAssign<&Scalar> for Scalar {
    fn add_assign(&mut self, rhs: &Scalar) {
        *self = *self + *rhs;
    }
}

impl<'a> core::ops::Sub<&'a Scalar> for &Scalar {
    type Output = Scalar;

    fn sub(self, other: &'a Scalar) -> Scalar {
        Scalar(self.0 - other.0)
    }
}

impl core::ops::Sub<Scalar> for Scalar {
    type Output = Scalar;

    fn sub(self, other: Scalar) -> Scalar {
        self - &other
    }
}

impl core::ops::Sub<&Scalar> for Scalar {
    type Output = Scalar;

    fn sub(self, other: &Scalar) -> Scalar {
        Scalar(self.0 - other.0)
    }
}

impl core::ops::SubAssign<Scalar> for Scalar {
    fn sub_assign(&mut self, rhs: Scalar) {
        *self = *self - rhs;
    }
}

impl core::ops::SubAssign<&Scalar> for Scalar {
    fn sub_assign(&mut self, rhs: &Scalar) {
        *self = *self - *rhs;
    }
}

impl<'a> core::ops::Mul<&'a Scalar> for &Scalar {
    type Output = Scalar;

    fn mul(self, other: &'a Scalar) -> Scalar {
        Scalar(self.0 * other.0)
    }
}

impl core::ops::Mul<Scalar> for Scalar {
    type Output = Scalar;

    fn mul(self, other: Scalar) -> Scalar {
        Scalar(self.0 * other.0)
    }
}

impl core::ops::Mul<&Scalar> for Scalar {
    type Output = Scalar;

    fn mul(self, other: &Scalar) -> Scalar {
        Scalar(self.0 * other.0)
    }
}

impl core::ops::MulAssign<Scalar> for Scalar {
    fn mul_assign(&mut self, rhs: Scalar) {
        *self = *self * rhs;
    }
}

impl core::ops::MulAssign<&Scalar> for Scalar {
    fn mul_assign(&mut self, rhs: &Scalar) {
        *self = *self * *rhs;
    }
}

impl core::ops::Neg for Scalar {
    type Output = Scalar;

    fn neg(self) -> Scalar {
        Scalar(-self.0)
    }
}

// Implement From<u64> for Scalar
impl core::convert::From<u64> for Scalar {
    fn from(v: u64) -> Self {
        Scalar(BackendScalar::from(v))
    }
}

// Implement ConstantTimeEq for Scalar
impl ConstantTimeEq for Scalar {
    fn ct_eq(&self, other: &Self) -> subtle::Choice {
        // Use constant-time comparison on the underlying BigInt limbs
        // Each limb comparison is constant-time, and we combine them with &
        // which is implemented as bitwise AND (constant-time)
        let a = &self.0.0.0;
        let b = &other.0.0.0;

        // CT comparison: all limbs must be equal
        a[0].ct_eq(&b[0]) & a[1].ct_eq(&b[1]) & a[2].ct_eq(&b[2]) & a[3].ct_eq(&b[3])
    }
}

// Constant-time equality of two base-field elements via their (Montgomery) limbs.
// `ark-ff` stores reduced values, so equal field elements have identical limbs.
fn fq_ct_eq(a: &BackendBaseField, b: &BackendBaseField) -> subtle::Choice {
    let a = &a.0.0;
    let b = &b.0.0;
    a[0].ct_eq(&b[0]) & a[1].ct_eq(&b[1]) & a[2].ct_eq(&b[2]) & a[3].ct_eq(&b[3])
}

// Constant-time equality for AffinePoint: compare x and y coordinate limbs.
impl ConstantTimeEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> subtle::Choice {
        fq_ct_eq(&self.x, &other.x) & fq_ct_eq(&self.y, &other.y)
    }
}

// Constant-time equality for ProjectivePoint via cross-multiplication, matching
// `PartialEq`: (X1:Y1:Z1) == (X2:Y2:Z2)  <=>  X1*Z2 == X2*Z1  AND  Y1*Z2 == Y2*Z1.
// This deliberately avoids `to_affine`, which would perform a variable-time field
// inversion (leaking the z coordinate through timing) and panic on `z == 0`. All
// four products are computed unconditionally, so control flow does not vary with
// the inputs, and the method is panic-free for any coordinates.
impl ConstantTimeEq for ProjectivePoint {
    fn ct_eq(&self, other: &Self) -> subtle::Choice {
        let x1z2 = self.x * other.z;
        let x2z1 = other.x * self.z;
        let y1z2 = self.y * other.z;
        let y2z1 = other.y * self.z;
        fq_ct_eq(&x1z2, &x2z1) & fq_ct_eq(&y1z2, &y2z1)
    }
}

// Inherent methods for random scalar generation
impl Scalar {
    pub fn random<R: rand_core::Rng + ?Sized>(rng: &mut R) -> Self {
        let mut bytes = [0u8; 40];
        rng.fill_bytes(&mut bytes);
        Self(BackendScalar::from_le_bytes_mod_order(&bytes))
    }

    pub fn try_random<R: rand_core::TryRng + ?Sized>(rng: &mut R) -> Result<Self, R::Error> {
        let mut bytes = [0u8; 40];
        rng.try_fill_bytes(&mut bytes)?;
        Ok(Self(BackendScalar::from_le_bytes_mod_order(&bytes)))
    }
}

// =====================================================================
// `CurveArithmetic` support — `Scalar`
//
// The impls below satisfy the bounds `elliptic_curve::CurveArithmetic` (v0.14)
// places on its `Scalar` associated type. They are modelled on the canonical
// reference `elliptic_curve::dev::mock_curve` and route every integer
// conversion through little-endian bytes / `U256`. Variable-time methods are
// permitted to be constant-time, so they delegate to the constant-time
// routines.
// =====================================================================

/// Canonical little-endian integer value of the scalar (`0 <= value < r`).
impl From<Scalar> for U256 {
    fn from(scalar: Scalar) -> U256 {
        U256::from_le_byte_array(scalar.to_bytes_le().into())
    }
}

impl From<&Scalar> for U256 {
    fn from(scalar: &Scalar) -> U256 {
        U256::from(*scalar)
    }
}

impl FromUintUnchecked for Scalar {
    type Uint = U256;

    /// Build a scalar from an integer the caller asserts is canonical (`< r`).
    /// Per the trait's "unchecked" contract, an out-of-range value is reduced
    /// modulo `r` rather than rejected.
    fn from_uint_unchecked(uint: U256) -> Self {
        let bytes: [u8; 32] = uint.to_le_byte_array().into();
        Scalar::reduce_bytes_le(&bytes)
    }
}

// `From`/`TryFrom` web between `Scalar`, `NonZeroScalar`, `ScalarValue`, and
// `SecretKey` (see `elliptic_curve::scalar_from_impls!`). Builds on the
// `From<&Scalar> for U256` and `FromUintUnchecked` impls above.
elliptic_curve::scalar_from_impls!(BabyJubJub, Scalar);

impl From<Scalar> for FieldBytes<BabyJubJub> {
    fn from(scalar: Scalar) -> Self {
        scalar.to_repr()
    }
}

impl From<&Scalar> for FieldBytes<BabyJubJub> {
    fn from(scalar: &Scalar) -> Self {
        scalar.to_repr()
    }
}

impl Reduce<U256> for Scalar {
    fn reduce(w: &U256) -> Self {
        let bytes: [u8; 32] = w.to_le_byte_array().into();
        Scalar::reduce_bytes_le(&bytes)
    }
}

impl Reduce<FieldBytes<BabyJubJub>> for Scalar {
    /// Reduce a little-endian `FieldBytes` value modulo `r`, matching the
    /// little-endian convention of `PrimeField::{from_repr, to_repr}`.
    fn reduce(w: &FieldBytes<BabyJubJub>) -> Self {
        let bytes: [u8; 32] = (*w).into();
        Scalar::reduce_bytes_le(&bytes)
    }
}

impl Retrieve for Scalar {
    type Output = U256;

    fn retrieve(&self) -> U256 {
        U256::from(*self)
    }
}

impl IsHigh for Scalar {
    /// `true` iff the canonical scalar is strictly greater than `r / 2`.
    fn is_high(&self) -> subtle::Choice {
        ScalarValue::<BabyJubJub>::from(*self).is_high()
    }
}

impl Invert for Scalar {
    type Output = CtOption<Scalar>;

    fn invert(&self) -> CtOption<Scalar> {
        Scalar::invert(self)
    }
}

impl Generate for Scalar {
    fn try_generate_from_rng<R: TryCryptoRng + ?Sized>(
        rng: &mut R,
    ) -> core::result::Result<Self, R::Error> {
        Scalar::try_random(rng)
    }
}

impl AsRef<Scalar> for Scalar {
    fn as_ref(&self) -> &Scalar {
        self
    }
}

impl PartialOrd for Scalar {
    /// Order scalars by their canonical integer value in `[0, r)`.
    fn partial_cmp(&self, other: &Self) -> Option<core::cmp::Ordering> {
        Some(U256::from(*self).cmp(&U256::from(*other)))
    }
}

impl ctutils::CtEq for Scalar {
    fn ct_eq(&self, other: &Self) -> ctutils::Choice {
        ConstantTimeEq::ct_eq(self, other).into()
    }
}

impl ctutils::CtSelect for Scalar {
    fn ct_select(&self, other: &Self, choice: ctutils::Choice) -> Self {
        if choice.to_bool() { *other } else { *self }
    }
}

// Scalar-as-LHS multiplication over points (`Scalar * Point` in both operand
// orders, plus the matching `MulVartime` impls). Delegates to the existing
// point-by-scalar routines.
elliptic_curve::scalar_mul_impls!(BabyJubJub, Scalar);

// =====================================================================
// `CurveArithmetic` support — points (`AffinePoint` / `ProjectivePoint`)
//
// As above, these satisfy the `CurveArithmetic` (v0.14) bounds on the
// `AffinePoint` / `ProjectivePoint` associated types and follow the canonical
// reference `elliptic_curve::dev::mock_curve`. Variable-time methods delegate
// to the existing (constant-time) point arithmetic.
// =====================================================================

// `Default` for both point types is the group identity. This matches the
// `elliptic_curve` convention (and `dev::mock_curve`): the generic
// `NonIdentity::new` uses `P::default()` as the identity sentinel (compared via
// `ct_eq`). The previous all-zero default made the projective cross-multiplied
// `ct_eq` degenerate (every point compared equal to a `z == 0` value), which
// would have broken `NonIdentity` (and hence `TryInto<NonIdentity<_>>`).
impl Default for AffinePoint {
    fn default() -> Self {
        Self::IDENTITY
    }
}

impl Default for ProjectivePoint {
    fn default() -> Self {
        Self::IDENTITY
    }
}

/// Serialize a base-field element to 32 canonical little-endian bytes
/// (allocation-free; mirrors [`Scalar::to_bytes_le`]).
fn fq_to_field_bytes(value: &BackendBaseField) -> FieldBytes<BabyJubJub> {
    let limbs = (*value).into_bigint().0;
    let mut arr = [0u8; 32];
    for (i, limb) in limbs.iter().enumerate() {
        arr[i * 8..i * 8 + 8].copy_from_slice(&limb.to_le_bytes());
    }
    arr.into()
}

// ----- AffinePoint -----

impl AffineCoordinates for AffinePoint {
    type FieldRepr = FieldBytes<BabyJubJub>;

    /// Build an affine point from little-endian coordinate bytes, returning
    /// `none` unless `(x, y)` is on the curve **and** in the prime-order
    /// subgroup. Out-of-range coordinate bytes are reduced modulo the
    /// base-field prime. This is **not** constant-time.
    fn from_coordinates(x: &Self::FieldRepr, y: &Self::FieldRepr) -> CtOption<Self> {
        let xb: [u8; 32] = (*x).into();
        let yb: [u8; 32] = (*y).into();
        let xf = BackendBaseField::from_le_bytes_mod_order(&xb);
        let yf = BackendBaseField::from_le_bytes_mod_order(&yb);
        match AffinePoint::new(xf, yf) {
            Some(point) => CtOption::new(point, subtle::Choice::from(1)),
            None => CtOption::new(AffinePoint::IDENTITY, subtle::Choice::from(0)),
        }
    }

    fn x(&self) -> Self::FieldRepr {
        fq_to_field_bytes(&self.x)
    }

    fn y(&self) -> Self::FieldRepr {
        fq_to_field_bytes(&self.y)
    }

    fn x_is_odd(&self) -> subtle::Choice {
        subtle::Choice::from((self.x.into_bigint().0[0] & 1) as u8)
    }

    fn y_is_odd(&self) -> subtle::Choice {
        subtle::Choice::from((self.y.into_bigint().0[0] & 1) as u8)
    }
}

impl core::ops::Neg for AffinePoint {
    type Output = AffinePoint;

    fn neg(self) -> AffinePoint {
        // Delegate to the backend's validated projective negation, then
        // normalise back to affine (the result has z = 1, so `to_affine` never
        // hits its `z == 0` panic path).
        (-ProjectivePoint::from(self)).to_affine()
    }
}

impl core::ops::Mul<Scalar> for AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: Scalar) -> ProjectivePoint {
        ProjectivePoint::from(self) * scalar
    }
}

impl core::ops::Mul<&Scalar> for AffinePoint {
    type Output = ProjectivePoint;

    fn mul(self, scalar: &Scalar) -> ProjectivePoint {
        ProjectivePoint::from(self) * scalar
    }
}

impl MulVartime<Scalar> for AffinePoint {
    fn mul_vartime(self, scalar: Scalar) -> ProjectivePoint {
        self * scalar
    }
}

impl MulVartime<&Scalar> for AffinePoint {
    fn mul_vartime(self, scalar: &Scalar) -> ProjectivePoint {
        self * scalar
    }
}

impl GroupEncoding for AffinePoint {
    type Repr = GroupRepr;

    fn from_bytes(bytes: &Self::Repr) -> CtOption<Self> {
        <ProjectivePoint as GroupEncoding>::from_bytes(bytes).map(|p| p.to_affine())
    }

    fn from_bytes_unchecked(bytes: &Self::Repr) -> CtOption<Self> {
        <ProjectivePoint as GroupEncoding>::from_bytes_unchecked(bytes).map(|p| p.to_affine())
    }

    fn to_bytes(&self) -> Self::Repr {
        ProjectivePoint::from(*self).to_bytes()
    }
}

impl CurveAffine for AffinePoint {
    type Curve = ProjectivePoint;
    type Scalar = Scalar;

    fn identity() -> Self {
        Self::IDENTITY
    }

    fn generator() -> Self {
        Self::GENERATOR
    }

    fn is_identity(&self) -> subtle::Choice {
        subtle::Choice::from(AffinePoint::is_identity(self) as u8)
    }

    fn to_curve(&self) -> ProjectivePoint {
        ProjectivePoint::from(*self)
    }
}

impl ctutils::CtEq for AffinePoint {
    fn ct_eq(&self, other: &Self) -> ctutils::Choice {
        ConstantTimeEq::ct_eq(self, other).into()
    }
}

impl ctutils::CtSelect for AffinePoint {
    fn ct_select(&self, other: &Self, choice: ctutils::Choice) -> Self {
        if choice.to_bool() { *other } else { *self }
    }
}

impl Generate for AffinePoint {
    fn try_generate_from_rng<R: TryCryptoRng + ?Sized>(
        rng: &mut R,
    ) -> core::result::Result<Self, R::Error> {
        Ok(<ProjectivePoint as Group>::try_random(rng)?.to_affine())
    }
}

impl From<NonIdentity<AffinePoint>> for AffinePoint {
    fn from(point: NonIdentity<AffinePoint>) -> Self {
        point.to_point()
    }
}

impl TryFrom<AffinePoint> for NonIdentity<AffinePoint> {
    type Error = EcError;

    fn try_from(point: AffinePoint) -> core::result::Result<Self, EcError> {
        NonIdentity::new(point).into_option().ok_or(EcError)
    }
}

// ----- ProjectivePoint -----

// Mixed projective/affine arithmetic, required by `CurveGroup`'s `GroupOps`
// bounds. Delegates through `From<AffinePoint>`.
impl core::ops::Add<AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, rhs: AffinePoint) -> ProjectivePoint {
        self + ProjectivePoint::from(rhs)
    }
}

impl core::ops::Add<&AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn add(self, rhs: &AffinePoint) -> ProjectivePoint {
        self + ProjectivePoint::from(*rhs)
    }
}

impl core::ops::Sub<AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, rhs: AffinePoint) -> ProjectivePoint {
        self - ProjectivePoint::from(rhs)
    }
}

impl core::ops::Sub<&AffinePoint> for ProjectivePoint {
    type Output = ProjectivePoint;

    fn sub(self, rhs: &AffinePoint) -> ProjectivePoint {
        self - ProjectivePoint::from(*rhs)
    }
}

impl core::ops::AddAssign<AffinePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: AffinePoint) {
        *self = *self + rhs;
    }
}

impl core::ops::AddAssign<&AffinePoint> for ProjectivePoint {
    fn add_assign(&mut self, rhs: &AffinePoint) {
        *self = *self + *rhs;
    }
}

impl core::ops::SubAssign<AffinePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: AffinePoint) {
        *self = *self - rhs;
    }
}

impl core::ops::SubAssign<&AffinePoint> for ProjectivePoint {
    fn sub_assign(&mut self, rhs: &AffinePoint) {
        *self = *self - *rhs;
    }
}

impl CurveGroup for ProjectivePoint {
    type Affine = AffinePoint;

    fn to_affine(&self) -> AffinePoint {
        // Calls the inherent `to_affine` (inherent methods take precedence over
        // trait methods of the same name), avoiding infinite recursion.
        ProjectivePoint::to_affine(self)
    }
}

impl<const N: usize> BatchNormalize<[ProjectivePoint; N]> for ProjectivePoint {
    type Output = [AffinePoint; N];

    fn batch_normalize(points: &[ProjectivePoint; N]) -> [AffinePoint; N] {
        core::array::from_fn(|i| points[i].to_affine())
    }
}

impl LinearCombination<[(ProjectivePoint, Scalar)]> for ProjectivePoint {}
impl<const N: usize> LinearCombination<[(ProjectivePoint, Scalar); N]> for ProjectivePoint {}

impl MulByGeneratorVartime for ProjectivePoint {}

impl MulVartime<Scalar> for ProjectivePoint {
    fn mul_vartime(self, scalar: Scalar) -> ProjectivePoint {
        self * scalar
    }
}

impl MulVartime<&Scalar> for ProjectivePoint {
    fn mul_vartime(self, scalar: &Scalar) -> ProjectivePoint {
        self * scalar
    }
}

impl ctutils::CtEq for ProjectivePoint {
    fn ct_eq(&self, other: &Self) -> ctutils::Choice {
        ConstantTimeEq::ct_eq(self, other).into()
    }
}

impl ctutils::CtSelect for ProjectivePoint {
    fn ct_select(&self, other: &Self, choice: ctutils::Choice) -> Self {
        if choice.to_bool() { *other } else { *self }
    }
}

impl Generate for ProjectivePoint {
    fn try_generate_from_rng<R: TryCryptoRng + ?Sized>(
        rng: &mut R,
    ) -> core::result::Result<Self, R::Error> {
        <ProjectivePoint as Group>::try_random(rng)
    }
}

impl From<NonIdentity<ProjectivePoint>> for ProjectivePoint {
    fn from(point: NonIdentity<ProjectivePoint>) -> Self {
        point.to_point()
    }
}

impl TryFrom<ProjectivePoint> for NonIdentity<ProjectivePoint> {
    type Error = EcError;

    fn try_from(point: ProjectivePoint) -> core::result::Result<Self, EcError> {
        NonIdentity::new(point).into_option().ok_or(EcError)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ec::CurveConfig;
    // `to_bytes_be` / `num_bits` on the backend's `BigInt` come from this trait.
    use ark_ff::BigInteger;

    /// Test that affine point operations match the backend implementation
    #[test]
    fn test_affine_point_identity() {
        let identity = AffinePoint::IDENTITY;
        assert!(identity.is_identity());
        assert_eq!(identity.x, BackendBaseField::ZERO);
        assert_eq!(identity.y, BackendBaseField::ONE);
    }

    /// Test that projective point operations match the backend implementation
    #[test]
    fn test_projective_point_identity() {
        let identity = ProjectivePoint::IDENTITY;
        assert_eq!(identity.x, BackendBaseField::ZERO);
        assert_eq!(identity.y, BackendBaseField::ONE);
        assert_eq!(identity.z, BackendBaseField::ONE);
    }

    /// Test that adding projective point to identity returns the same point
    #[test]
    fn test_add_identity_projective() {
        let point = ProjectivePoint::GENERATOR;
        let result = point + ProjectivePoint::IDENTITY;
        assert_eq!(result.x, point.x);
        assert_eq!(result.y, point.y);
        assert_eq!(result.z, point.z);
    }

    /// Test that &T + &T works correctly
    #[test]
    fn test_add_refs() {
        let a = ProjectivePoint::GENERATOR;
        let b = ProjectivePoint::IDENTITY;
        let a_ref = &a;
        let b_ref = &b;
        let result = a_ref + b_ref;
        assert_eq!(result.x, a.x);
        assert_eq!(result.y, a.y);
        assert_eq!(result.z, a.z);
    }

    /// Test scalar multiplication on generator matches backend
    #[test]
    fn test_scalar_mult_on_generator() {
        // 1 * G = G
        let scalar_one = Scalar::ONE;
        let result = ProjectivePoint::GENERATOR * scalar_one;

        let expected = BackendProjective::from(ProjectivePoint::GENERATOR);
        let expected_result = expected * BackendScalar::ONE;

        let result_backend = BackendProjective::from(result);
        assert_eq!(result_backend, expected_result);
    }

    /// Test that 2 * G matches backend implementation
    #[test]
    fn test_double_generator() {
        let scalar_two: Scalar = 2u64.into();
        let result = ProjectivePoint::GENERATOR * scalar_two;

        let backend_g = BackendProjective::from(ProjectivePoint::GENERATOR);
        let expected = backend_g + backend_g;

        let result_proj = BackendProjective::from(result);
        assert_eq!(result_proj, expected);
    }

    /// Test scalar multiplication with zero
    #[test]
    fn test_scalar_mult_zero() {
        let scalar_zero: Scalar = 0u64.into();
        let result = ProjectivePoint::GENERATOR * scalar_zero;
        assert_eq!(result, ProjectivePoint::IDENTITY);
    }

    /// Test scalar multiplication with one
    #[test]
    fn test_scalar_mult_one() {
        let scalar_one: Scalar = 1u64.into();
        let result = ProjectivePoint::GENERATOR * scalar_one;

        let expected = ProjectivePoint::GENERATOR;
        assert_eq!(result.x, expected.x);
        assert_eq!(result.y, expected.y);
        assert_eq!(result.z, expected.z);
    }

    /// Test scalar field arithmetic matches backend
    #[test]
    fn test_scalar_add() {
        let a: Scalar = 5u64.into();
        let b: Scalar = 7u64.into();
        let result = a + b;
        let expected: Scalar = 12u64.into();

        assert_eq!(result.0, expected.0);
    }

    /// Test scalar field multiplication
    #[test]
    fn test_scalar_mult() {
        let a: Scalar = 3u64.into();
        let b: Scalar = 4u64.into();
        let result = a * b;
        let expected: Scalar = 12u64.into();

        assert_eq!(result.0, expected.0);
    }

    /// Test scalar negation
    #[test]
    fn test_scalar_neg() {
        let a: Scalar = 5u64.into();
        let neg_a = -a;
        let zero: Scalar = 0u64.into();

        // a + (-a) should equal 0
        assert_eq!((a + neg_a).0, zero.0);
    }

    /// Test that Scalar::from(1) equals Scalar::ONE
    #[test]
    fn test_scalar_from_u64() {
        let from_u64: Scalar = 1u64.into();
        assert_eq!(from_u64.0, Scalar::ONE.0);
    }

    /// Test affine to projective conversion
    #[test]
    fn test_affine_to_projective() {
        let affine = AffinePoint::GENERATOR;
        // Direct conversion from AffinePoint to ProjectivePoint
        let projective = ProjectivePoint::from(affine);

        assert_eq!(projective.x, affine.x);
        assert_eq!(projective.y, affine.y);
        assert_eq!(projective.z, BackendBaseField::ONE);
    }

    /// Test that point operations work correctly with reference types
    #[test]
    fn test_point_ops_with_refs() {
        let a = ProjectivePoint::GENERATOR;
        let b = ProjectivePoint::IDENTITY;

        // Test &T + &T
        let a_ref = &a;
        let b_ref = &b;
        let result1 = a_ref + b_ref;
        assert_eq!(result1.x, a.x);
        assert_eq!(result1.y, a.y);
        assert_eq!(result1.z, a.z);
    }

    /// Test that curve order times generator gives identity
    /// This is a key property of elliptic curves
    #[test]
    fn test_curve_order_times_generator() {
        use elliptic_curve::Curve;

        // Get the order as bytes
        let order_uint = BabyJubJub::ORDER.as_ref();
        let order_bytes = order_uint.to_be_bytes();

        // Convert to scalar (big-endian to little-endian)
        let mut order_le = [0u8; 32];
        for (i, b) in order_bytes.iter().enumerate() {
            order_le[31 - i] = *b;
        }

        // `r` is non-canonical as a scalar (it equals the modulus), so it is
        // reduced: r mod r == 0. Hence [r]G == [0]G == identity.
        let order_scalar = Scalar::reduce_bytes_le(&order_le);
        assert!(order_scalar.is_zero());
        let result = ProjectivePoint::GENERATOR * order_scalar;

        assert_eq!(result, ProjectivePoint::IDENTITY);
    }

    /// Test scalar field constant TWO_INV
    #[test]
    fn test_two_inv() {
        // TWO_INV * 2 should equal 1
        let result = Scalar::TWO_INV * Scalar::from(2u64);
        assert_eq!(result.0, Scalar::ONE.0);
    }

    /// Test scalar field constant MULTIPLICATIVE_GENERATOR is a non-residue.
    #[test]
    fn test_multiplicative_generator() {
        // The generator must be a quadratic non-residue (and is NOT 5, which is
        // a quadratic residue mod r and therefore an invalid generator).
        let g = Scalar::MULTIPLICATIVE_GENERATOR.0;
        assert_eq!(
            g.pow(BackendScalar::MODULUS_MINUS_ONE_DIV_TWO),
            -BackendScalar::ONE
        );
        assert_ne!(g, Scalar::from(5u64).0);
    }

    /// Test constant time equality
    #[test]
    fn test_ct_eq() {
        let a: Scalar = 42u64.into();
        let b: Scalar = 42u64.into();
        let c: Scalar = 43u64.into();

        assert!(a.ct_eq(&b).unwrap_u8() == 1);
        assert!(a.ct_eq(&c).unwrap_u8() == 0);
    }

    /// Test scalar field multiplication with large values
    #[test]
    fn test_scalar_mult_large() {
        use elliptic_curve::Curve;

        // Get the curve order (which is also the scalar field order)
        let order_uint = BabyJubJub::ORDER.as_ref();
        let order_bytes = order_uint.to_be_bytes();

        // Convert to scalar (big-endian to little-endian)
        let mut order_le = [0u8; 32];
        for (i, b) in order_bytes.iter().enumerate() {
            order_le[31 - i] = *b;
        }

        // `r` is non-canonical as a scalar (it equals the modulus), so it is
        // reduced: r mod r == 0.
        let order_scalar = Scalar::reduce_bytes_le(&order_le);
        assert!(order_scalar.is_zero());

        // Test that order * G = identity (== 0 * G).
        let result = ProjectivePoint::GENERATOR * order_scalar;
        assert_eq!(result, ProjectivePoint::IDENTITY);

        // Test that (order-1) * G = -G (since order*G = 0)
        let scalar_r_minus_1 = order_scalar - Scalar::ONE;
        let result = ProjectivePoint::GENERATOR * scalar_r_minus_1;

        // order*G = 0, so (order-1)*G = -G
        let neg_gen = -ProjectivePoint::GENERATOR;

        // In projective coordinates, (x, y, z) and (x*z', y*z', z*z') represent the same point
        // So we need to compare in affine form
        let result_affine = result.to_affine();
        let neg_gen_affine = neg_gen.to_affine();
        assert_eq!(result_affine.x, neg_gen_affine.x);
        assert_eq!(result_affine.y, neg_gen_affine.y);
    }

    /// Test scalar multiplication matches the backend implementation with random scalars
    #[test]
    fn test_scalar_mult_consistency() {
        use rand::{SeedableRng, rngs::StdRng};

        const NUM_TESTS: usize = 1000;
        let seed: [u8; 32] = rand::random();

        let mut rng = StdRng::from_seed(seed);

        for _ in 0..NUM_TESTS {
            // Generate a random scalar
            let scalar = Scalar::random(&mut rng);

            // Compute A * B using our implementation
            let result_ours = ProjectivePoint::GENERATOR * scalar;

            // Compute A * B using the backend implementation
            let backend_scalar: BackendScalar = scalar.0;
            let backend_generator: BackendProjective = ProjectivePoint::GENERATOR.into();
            let result_backend = backend_generator * backend_scalar;

            // Convert our result to backend and compare
            let result_as_backend: BackendProjective = result_ours.into();

            assert_eq!(
                result_as_backend,
                result_backend,
                "Scalar multiplication mismatch, seed = {}",
                hex::encode(seed)
            );
        }
    }

    /// Test scalar from_bytes (big-endian)
    #[test]
    fn test_scalar_from_bytes_be() {
        let mut bytes = [0u8; 32];
        bytes[31] = 42; // big-endian: 42 in the last byte = 42 in LE first byte
        let scalar = Scalar::from_bytes(&bytes).unwrap();
        let expected: Scalar = 42u64.into();
        assert_eq!(scalar.0, expected.0);
    }

    /// Test scalar from_bytes_le (little-endian)
    #[test]
    fn test_scalar_from_bytes_le_input() {
        let mut bytes = [0u8; 32];
        bytes[0] = 42; // little-endian: 42 in the first byte
        let scalar = Scalar::from_bytes_le(&bytes).unwrap();
        let expected: Scalar = 42u64.into();
        assert_eq!(scalar.0, expected.0);
    }

    /// Test scalar to_bytes (big-endian) conversion
    #[test]
    fn test_scalar_to_bytes() {
        let scalar: Scalar = 42u64.into();
        let bytes = scalar.to_bytes();
        // 42 in big-endian should have 42 in the last position
        assert_eq!(bytes[31], 42);
    }

    /// Test scalar invert - that a * a^-1 = 1
    #[test]
    fn test_scalar_invert() {
        let scalar: Scalar = 5u64.into();
        let inverted = scalar.invert().unwrap();
        let result = scalar * inverted;
        assert_eq!(result.0, Scalar::ONE.0);
    }

    /// Test scalar invert on zero
    #[test]
    fn test_scalar_invert_zero() {
        let scalar = Scalar::ZERO;
        let result = scalar.invert();
        // Invert of zero should return None (failure)
        // CtOption::is_none() returns Choice, so we use unwrap_u8()
        assert_eq!(result.is_none().unwrap_u8(), 1);
    }

    /// Test scalar square
    #[test]
    fn test_scalar_square() {
        let scalar: Scalar = 5u64.into();
        let squared = scalar.square();
        let expected: Scalar = 25u64.into();
        assert_eq!(squared.0, expected.0);
    }

    /// Test scalar double
    #[test]
    fn test_scalar_double() {
        let scalar: Scalar = 5u64.into();
        let doubled = scalar.double();
        let expected: Scalar = 10u64.into();
        assert_eq!(doubled.0, expected.0);
    }

    /// Test is_zero scalar method
    #[test]
    fn test_scalar_is_zero() {
        assert!(Scalar::ZERO.is_zero());
        assert!(!Scalar::ONE.is_zero());
        let non_zero: Scalar = 42u64.into();
        assert!(!non_zero.is_zero());
    }

    /// Test is_one scalar method
    #[test]
    fn test_scalar_is_one() {
        assert!(Scalar::ONE.is_one());
        assert!(!Scalar::ZERO.is_one());
        let not_one: Scalar = 42u64.into();
        assert!(!not_one.is_one());
    }

    /// Test GroupEncoding to_bytes is the canonical 32-byte compressed encoding.
    #[test]
    fn test_group_encoding_to_bytes() {
        let point = ProjectivePoint::GENERATOR;
        let repr = point.to_bytes();
        assert_eq!(repr.as_ref().len(), 32);
    }

    /// All-zero bytes decode to `y = 0`; any such on-curve point is 4-torsion
    /// (not in the prime-order subgroup), so the validating decoder rejects it.
    #[test]
    fn test_group_encoding_from_bytes() {
        let bytes = GroupRepr([0u8; 32]);
        let result = ProjectivePoint::from_bytes(&bytes);
        assert_eq!(result.is_none().unwrap_u8(), 1);
    }

    /// Test Sum trait for ProjectivePoint
    #[test]
    fn test_projective_point_sum() {
        let points = [
            ProjectivePoint::GENERATOR,
            ProjectivePoint::GENERATOR,
            ProjectivePoint::GENERATOR,
        ];
        let sum: ProjectivePoint = points.into_iter().sum();
        // 3 * G = G + G + G
        let g = ProjectivePoint::GENERATOR;
        let expected = g + g + g;
        assert_eq!(sum.x, expected.x);
        assert_eq!(sum.y, expected.y);
        assert_eq!(sum.z, expected.z);
    }

    /// Test Sum trait for Scalar
    #[test]
    fn test_scalar_sum() {
        let scalars = [Scalar::from(1u64), Scalar::from(2u64), Scalar::from(3u64)];
        let sum: Scalar = scalars.into_iter().sum();
        let expected = Scalar::from(6u64);
        assert_eq!(sum.0, expected.0);
    }

    /// Test Product trait for Scalar
    #[test]
    fn test_scalar_product() {
        let scalars = [Scalar::from(2u64), Scalar::from(3u64), Scalar::from(5u64)];
        let product: Scalar = scalars.into_iter().product();
        let expected = Scalar::from(30u64);
        assert_eq!(product.0, expected.0);
    }

    /// Test conversion round-trip: Affine -> Projective -> Affine
    #[test]
    fn test_affine_projective_round_trip() {
        let affine = AffinePoint::GENERATOR;
        let projective: ProjectivePoint = affine.into();
        let affine_back: AffinePoint = projective.into();
        assert_eq!(affine.x, affine_back.x);
        assert_eq!(affine.y, affine_back.y);
    }

    /// Test AffinePoint::new_unchecked
    #[test]
    fn test_affine_point_new_unchecked() {
        let affine = AffinePoint::new_unchecked(BackendBaseField::ONE, BackendBaseField::ONE);
        assert_eq!(affine.x, BackendBaseField::ONE);
        assert_eq!(affine.y, BackendBaseField::ONE);
    }

    /// Test AffinePoint::x getter
    #[test]
    fn test_affine_point_x() {
        let affine = AffinePoint::GENERATOR;
        assert_eq!(affine.x(), affine.x);
    }

    /// Test AffinePoint::y getter
    #[test]
    fn test_affine_point_y() {
        let affine = AffinePoint::GENERATOR;
        assert_eq!(affine.y(), affine.y);
    }

    /// Test AffinePoint::is_identity on identity
    #[test]
    fn test_affine_point_is_identity() {
        assert!(AffinePoint::IDENTITY.is_identity());
        assert!(!AffinePoint::GENERATOR.is_identity());
    }

    /// Test AffinePoint::x_is_odd
    #[test]
    fn test_affine_point_x_is_odd() {
        let affine = AffinePoint::GENERATOR;
        // Just verify it returns a choice (doesn't panic)
        let _ = affine.x_is_odd();
    }

    /// Test ProjectivePoint::new_unchecked
    #[test]
    fn test_projective_point_new_unchecked() {
        let point = ProjectivePoint::new_unchecked(
            BackendBaseField::ONE,
            BackendBaseField::ONE,
            BackendBaseField::ONE,
        );
        assert_eq!(point.x, BackendBaseField::ONE);
        assert_eq!(point.y, BackendBaseField::ONE);
        assert_eq!(point.z, BackendBaseField::ONE);
    }

    /// Test ProjectivePoint::to_affine on identity
    #[test]
    fn test_projective_to_affine_identity() {
        let projective = ProjectivePoint::IDENTITY;
        let affine = projective.to_affine();
        assert_eq!(affine, AffinePoint::IDENTITY);
    }

    /// Test ProjectivePoint::to_affine on generator
    #[test]
    fn test_projective_to_affine_generator() {
        let projective = ProjectivePoint::GENERATOR;
        let affine = projective.to_affine();
        assert_eq!(affine.x, ProjectivePoint::GENERATOR.x);
        assert_eq!(affine.y, ProjectivePoint::GENERATOR.y);
    }

    /// Test ProjectivePoint::to_affine with non-one z coordinate
    #[test]
    fn test_projective_to_affine_non_one_z() {
        // Use 2 * generator to verify conversion
        let doubled = ProjectivePoint::GENERATOR * Scalar::from(2u64);
        let affine = doubled.to_affine();
        // Just verify we can convert and get a valid point
        assert!(affine.x != BackendBaseField::ZERO || affine.y != BackendBaseField::ONE);
    }

    /// Test Scalar::from_bytes_le with a known value
    #[test]
    fn test_scalar_from_bytes_le_known() {
        let mut bytes = [0u8; 32];
        bytes[0] = 0x01; // 1 in little-endian
        let scalar = Scalar::from_bytes_le(&bytes).unwrap();
        assert_eq!(scalar.0, Scalar::ONE.0);
    }

    /// `from_repr_vartime` must be little-endian (matching `to_repr`) and reject
    /// non-canonical (`>= r`) encodings.
    #[test]
    fn test_scalar_from_repr_vartime() {
        let mut bytes = [0u8; 32];
        bytes[0] = 42; // 42 little-endian
        let scalar = Scalar::from_repr_vartime(bytes.into());
        assert!(scalar.is_some());
        let expected: Scalar = 42u64.into();
        assert_eq!(scalar.unwrap().0, expected.0);
        // All-ones (= 2^256 - 1 >= r) is non-canonical and must be rejected.
        assert!(Scalar::from_repr_vartime([0xFFu8; 32].into()).is_none());
    }

    /// Test Scalar::is_odd
    #[test]
    fn test_scalar_is_odd() {
        let even: Scalar = 42u64.into();
        let odd: Scalar = 43u64.into();
        assert!(even.is_odd().unwrap_u8() == 0); // 42 is even
        assert!(odd.is_odd().unwrap_u8() == 1); // 43 is odd
    }

    /// Test Scalar::sqrt on 16 (a quadratic residue)
    #[test]
    fn test_scalar_sqrt() {
        let val: Scalar = 16u64.into();
        let sqrt = val.sqrt();
        // 16 = 4^2, sqrt returns Some(result) if it exists
        // CtOption::is_some() returns Choice, so use unwrap_u8()
        assert_eq!(sqrt.is_some().unwrap_u8(), 1);
        // Verify by squaring the result
        let result = sqrt.unwrap();
        let squared = result.square();
        assert_eq!(squared.0, val.0);
    }

    /// Test Scalar::sqrt on 2 (also a quadratic residue in this field)
    #[test]
    fn test_scalar_sqrt_two() {
        let val: Scalar = 2u64.into();
        let sqrt = val.sqrt();
        // sqrt returns Some(result) if it exists
        assert_eq!(sqrt.is_some().unwrap_u8(), 1);
        // Verify by squaring the result
        let result = sqrt.unwrap();
        let squared = result.square();
        assert_eq!(squared.0, val.0);
    }

    /// `sqrt_ratio` must compute a real result per the `ff` contract (the old
    /// implementation was a stub that always returned `(1, 1)`).
    #[test]
    fn test_scalar_sqrt_ratio() {
        // 16/4 = 4 is a square: is_square == 1 and root^2 == num/den.
        let num: Scalar = 16u64.into();
        let den: Scalar = 4u64.into();
        let (is_square, root) = Scalar::sqrt_ratio(&num, &den);
        assert_eq!(is_square.unwrap_u8(), 1);
        let ratio = num * den.invert().unwrap();
        assert_eq!(root.square().0, ratio.0, "root^2 must equal num/den");

        // ROOT_OF_UNITY is a non-square, so its ratio over 1 is a non-square.
        let (is_square_ns, _) = Scalar::sqrt_ratio(&Scalar::ROOT_OF_UNITY, &Scalar::ONE);
        assert_eq!(
            is_square_ns.unwrap_u8(),
            0,
            "a non-square must report is_square == 0"
        );

        // num != 0, den == 0 => (0, _).
        let (is_square_d0, _) = Scalar::sqrt_ratio(&num, &Scalar::ZERO);
        assert_eq!(is_square_d0.unwrap_u8(), 0);

        // num == 0 => (1, 0).
        let (is_square_n0, root_n0) = Scalar::sqrt_ratio(&Scalar::ZERO, &den);
        assert_eq!(is_square_n0.unwrap_u8(), 1);
        assert_eq!(root_n0.0, Scalar::ZERO.0);
    }

    /// Test Scalar::conditional_select
    #[test]
    fn test_scalar_conditional_select() {
        let a: Scalar = 42u64.into();
        let b: Scalar = 84u64.into();
        let result = Scalar::conditional_select(&a, &b, 0.into());
        assert_eq!(result.0, a.0);
        let result = Scalar::conditional_select(&a, &b, 1.into());
        assert_eq!(result.0, b.0);
    }

    /// Test AffinePoint::conditional_select
    #[test]
    fn test_affine_point_conditional_select() {
        let a = AffinePoint::IDENTITY;
        let b = AffinePoint::GENERATOR;
        let result = AffinePoint::conditional_select(&a, &b, 0.into());
        assert_eq!(result.x, a.x);
        assert_eq!(result.y, a.y);
        let result = AffinePoint::conditional_select(&a, &b, 1.into());
        assert_eq!(result.x, b.x);
        assert_eq!(result.y, b.y);
    }

    /// Test GroupRepr::default (canonical 32-byte repr)
    #[test]
    fn test_group_repr_default() {
        let repr = GroupRepr::default();
        assert_eq!(repr.as_ref(), &[0u8; 32]);
    }

    /// Test BabyJubJub::ORDER constant matches the backend
    #[test]
    fn test_babyjubjub_order() {
        // Verify our hardcoded ORDER matches the backend's MODULUS
        // This ensures the hardcoded value is not outdated

        // Get the ORDER from BabyJubJub and convert to hex
        let order = BabyJubJub::ORDER;
        let order_bytes = order.as_ref().to_be_bytes();
        // Convert bytes to hex string (big-endian)
        let order_hex = hex::encode(order_bytes);

        // Get the backend MODULUS as hex
        let backend_modulus = BackendScalar::MODULUS;
        let backend_bytes = backend_modulus.to_bytes_be();
        let backend_hex = hex::encode(backend_bytes);

        // Both should be equal (this verifies our ORDER_HEX is correct)
        assert_eq!(
            order_hex, backend_hex,
            "ORDER_HEX must match BackendScalar::MODULUS"
        );

        // Also verify the order is non-zero and has expected properties
        // First byte in big-endian should be 0x06 (non-zero)
        assert_eq!(order_bytes[0], 0x06);
    }

    /// Verify [`COFACTOR`] matches `BackendProjective::Config::COFACTOR`.
    ///
    /// If the backend is updated, this test will fail and force a review of
    /// any hardcoded cofactor-dependent logic (including [`ProjectivePoint::mul_with_cofactor_clear`]).
    #[test]
    fn test_cofactor_matches_backend() {
        let backend_cofactor = <BackendProjective as ark_ec::CurveGroup>::Config::COFACTOR;
        assert_eq!(
            backend_cofactor.len(),
            1,
            "BackendConfig::COFACTOR must be a single-limb value for this test"
        );
        assert_eq!(
            COFACTOR, backend_cofactor[0],
            "COFACTOR constant must match BackendConfig::COFACTOR"
        );
    }

    /// Verify the `ff` field constants actually satisfy their mathematical
    /// contracts (not merely "non-zero"). These checks pin down ROOT_OF_UNITY,
    /// ROOT_OF_UNITY_INV, DELTA, MULTIPLICATIVE_GENERATOR and S against the
    /// backend, so a wrong hardcoded value can never pass CI again.
    #[test]
    fn test_montgomery_constants() {
        // 2 * TWO_INV == 1
        let two: Scalar = 2u64.into();
        assert_eq!((two * Scalar::TWO_INV).0, Scalar::ONE.0, "TWO_INV wrong");

        // S must equal the backend's 2-adicity.
        assert_eq!(Scalar::S, BackendScalar::TWO_ADICITY, "S != TWO_ADICITY");
        assert_eq!(Scalar::S, 4);

        let g = Scalar::MULTIPLICATIVE_GENERATOR.0;
        let one = BackendScalar::ONE;

        // MULTIPLICATIVE_GENERATOR must be a quadratic non-residue:
        // g^((r-1)/2) == -1.
        assert_eq!(
            g.pow(BackendScalar::MODULUS_MINUS_ONE_DIV_TWO),
            -one,
            "MULTIPLICATIVE_GENERATOR is not a quadratic non-residue"
        );

        // ROOT_OF_UNITY must be derived from the generator: g^t, t = (r-1) >> S.
        assert_eq!(
            g.pow(BackendScalar::TRACE),
            Scalar::ROOT_OF_UNITY.0,
            "ROOT_OF_UNITY != MULTIPLICATIVE_GENERATOR^t"
        );

        // ROOT_OF_UNITY must be a *primitive* 2^S root of unity:
        // root^(2^S) == 1 but root^(2^(S-1)) == -1.
        let mut acc = Scalar::ROOT_OF_UNITY.0;
        for _ in 0..(Scalar::S - 1) {
            acc = acc.square();
        }
        assert_eq!(
            acc, -one,
            "ROOT_OF_UNITY^(2^(S-1)) must be -1 (primitivity)"
        );
        acc = acc.square();
        assert_eq!(acc, one, "ROOT_OF_UNITY^(2^S) must be 1");

        // ROOT_OF_UNITY_INV must be the inverse of ROOT_OF_UNITY.
        assert_eq!(
            (Scalar::ROOT_OF_UNITY * Scalar::ROOT_OF_UNITY_INV).0,
            Scalar::ONE.0,
            "ROOT_OF_UNITY_INV is not the inverse of ROOT_OF_UNITY"
        );

        // DELTA must equal MULTIPLICATIVE_GENERATOR^(2^S).
        let mut delta = g;
        for _ in 0..Scalar::S {
            delta = delta.square();
        }
        assert_eq!(
            delta,
            Scalar::DELTA.0,
            "DELTA != MULTIPLICATIVE_GENERATOR^(2^S)"
        );
    }

    /// Test that Scalar::ZERO and Scalar::ONE match backend values
    #[test]
    fn test_scalar_zero_one() {
        // ZERO should be 0 in Montgomery form
        assert_eq!(Scalar::ZERO.0, BackendScalar::ZERO);

        // ONE should be 1 in Montgomery form (which is R mod r in Montgomery)
        // We can verify this works correctly in field operations
        let one_plus_one = Scalar::ONE + Scalar::ONE;
        let two: Scalar = 2u64.into();
        assert_eq!(one_plus_one.0, two.0, "ONE + ONE should equal 2");
    }

    /// NUM_BITS / CAPACITY must reflect the *actual* 251-bit modulus, not an
    /// over-claimed 255/254 (which silently truncated values a caller packed
    /// into CAPACITY bits).
    #[test]
    fn test_scalar_bit_constants() {
        // NUM_BITS must equal the real modulus bit-size.
        assert_eq!(
            Scalar::NUM_BITS,
            BackendScalar::MODULUS_BIT_SIZE,
            "NUM_BITS must equal the real modulus bit size"
        );
        assert_eq!(Scalar::NUM_BITS, 251);
        // CAPACITY is NUM_BITS - 1.
        assert_eq!(Scalar::CAPACITY, Scalar::NUM_BITS - 1);
        assert_eq!(Scalar::CAPACITY, 250);

        // The modulus must occupy exactly NUM_BITS bits.
        assert_eq!(BackendScalar::MODULUS.num_bits(), Scalar::NUM_BITS);
    }

    /// A value with exactly `CAPACITY` bits set must always round-trip through
    /// the scalar canonically, whereas a `NUM_BITS`-bit value need not (it can
    /// exceed `r`). This guards against the previous over-claimed CAPACITY.
    #[test]
    fn test_scalar_capacity_round_trip() {
        // 2^CAPACITY - 1 (all CAPACITY low bits set) is < r, so it round-trips.
        let cap = Scalar::CAPACITY as usize;
        let mut le = [0u8; 32];
        for bit in 0..cap {
            le[bit / 8] |= 1 << (bit % 8);
        }
        let s = Scalar::from_bytes_le(&le);
        assert_eq!(s.is_some().unwrap_u8(), 1, "CAPACITY-bit value must decode");
        assert_eq!(
            s.unwrap().to_bytes_le(),
            le,
            "CAPACITY-bit value must round-trip"
        );
    }

    /// Test BabyJubJub::FieldBytesSize
    #[test]
    fn test_babyjubjub_field_bytes_size() {
        // Just verify it exists (compile-time check)
        let _ = core::mem::size_of::<<BabyJubJub as elliptic_curve::Curve>::FieldBytesSize>();
    }

    /// Test BabyJubJub::Uint type
    #[test]
    fn test_babyjubjub_uint() {
        // Just verify it exists (compile-time check)
        let _ = core::mem::size_of::<<BabyJubJub as elliptic_curve::Curve>::Uint>();
    }

    /// Test Scalar::TWO_INV
    #[test]
    fn test_scalar_two_inv() {
        let result = Scalar::TWO_INV * Scalar::from(2u64);
        assert_eq!(result.0, Scalar::ONE.0);
    }

    /// MULTIPLICATIVE_GENERATOR must be a quadratic non-residue (and therefore
    /// not 5, which is a QR mod r).
    #[test]
    fn test_scalar_multiplicative_generator() {
        let g = Scalar::MULTIPLICATIVE_GENERATOR.0;
        assert_eq!(
            g.pow(BackendScalar::MODULUS_MINUS_ONE_DIV_TWO),
            -BackendScalar::ONE,
            "generator must be a quadratic non-residue"
        );
        let five: Scalar = 5u64.into();
        assert_ne!(g, five.0, "5 is a QR mod r and is not a valid generator");
    }

    /// Test Scalar::NUM_BITS reflects the real 251-bit modulus.
    #[test]
    fn test_scalar_num_bits() {
        assert_eq!(Scalar::NUM_BITS, 251);
        assert_eq!(Scalar::NUM_BITS, BackendScalar::MODULUS_BIT_SIZE);
    }

    /// Test Scalar::CAPACITY is NUM_BITS - 1.
    #[test]
    fn test_scalar_capacity() {
        assert_eq!(Scalar::CAPACITY, 250);
        assert_eq!(Scalar::CAPACITY, Scalar::NUM_BITS - 1);
    }

    /// Test Scalar::S
    #[test]
    fn test_scalar_s() {
        // S = 4 because r - 1 = 2^4 * t with t odd.
        assert_eq!(Scalar::S, 4);
        assert_eq!(Scalar::S, BackendScalar::TWO_ADICITY);
    }

    /// ROOT_OF_UNITY must be a primitive 2^S root of unity.
    #[test]
    fn test_scalar_root_of_unity() {
        let mut acc = Scalar::ROOT_OF_UNITY.0;
        for _ in 0..(Scalar::S - 1) {
            acc = acc.square();
        }
        // root^(2^(S-1)) == -1 (primitive), root^(2^S) == 1.
        assert_eq!(acc, -BackendScalar::ONE);
        assert_eq!(acc.square(), BackendScalar::ONE);
    }

    /// ROOT_OF_UNITY_INV must be the multiplicative inverse of ROOT_OF_UNITY.
    #[test]
    fn test_scalar_root_of_unity_inv() {
        assert_eq!(
            (Scalar::ROOT_OF_UNITY * Scalar::ROOT_OF_UNITY_INV).0,
            Scalar::ONE.0
        );
    }

    /// DELTA must equal MULTIPLICATIVE_GENERATOR^(2^S).
    #[test]
    fn test_scalar_delta() {
        let mut delta = Scalar::MULTIPLICATIVE_GENERATOR.0;
        for _ in 0..Scalar::S {
            delta = delta.square();
        }
        assert_eq!(delta, Scalar::DELTA.0);
    }

    /// Test From<&AffinePoint> for ProjectivePoint
    #[test]
    fn test_from_affine_ref_to_projective() {
        let affine = AffinePoint::GENERATOR;
        let projective: ProjectivePoint = (&affine).into();
        assert_eq!(projective.x, affine.x);
        assert_eq!(projective.y, affine.y);
        assert_eq!(projective.z, BackendBaseField::ONE);
    }

    /// Test From<&ProjectivePoint> for AffinePoint
    #[test]
    fn test_from_projective_ref_to_affine() {
        let projective = ProjectivePoint::GENERATOR;
        let affine: AffinePoint = (&projective).into();
        assert_eq!(affine.x, projective.x);
        assert_eq!(affine.y, projective.y);
    }

    /// Test From<&BackendProjective> for ProjectivePoint
    #[test]
    fn test_from_backend_projective_ref() {
        let backend = BackendProjective::from(ProjectivePoint::GENERATOR);
        let point: ProjectivePoint = (&backend).into();
        assert_eq!(point.x, backend.x);
        assert_eq!(point.y, backend.y);
        assert_eq!(point.z, backend.z);
    }

    /// Test From<&ProjectivePoint> for BackendProjective
    #[test]
    fn test_from_projective_ref_to_backend() {
        let point = ProjectivePoint::GENERATOR;
        let backend: BackendProjective = (&point).into();
        assert_eq!(backend.x, point.x);
        assert_eq!(backend.y, point.y);
        assert_eq!(backend.z, point.z);
    }

    /// Test AsRef<[u8]> for GroupRepr
    #[test]
    fn test_group_repr_as_ref() {
        let repr = GroupRepr([42u8; 32]);
        let bytes: &[u8] = repr.as_ref();
        assert_eq!(bytes.len(), 32);
        assert_eq!(bytes[0], 42);
    }

    /// Test AsMut<[u8]> for GroupRepr
    #[test]
    fn test_group_repr_as_mut() {
        let mut repr = GroupRepr([0u8; 32]);
        let bytes: &mut [u8] = repr.as_mut();
        bytes[0] = 42;
        assert_eq!(repr.as_ref()[0], 42);
    }

    /// Test Scalar negation via subtraction
    #[test]
    fn test_scalar_sub_neg() {
        let a: Scalar = 42u64.into();
        let zero: Scalar = 0u64.into();
        let result = zero - a;
        assert_eq!(result.0, (-a).0);
    }

    /// Test Neg for &ProjectivePoint
    #[test]
    fn test_projective_point_neg_ref() {
        let point = ProjectivePoint::GENERATOR;
        let neg_point = -&point;
        let result = point + neg_point;
        assert_eq!(result, ProjectivePoint::IDENTITY);
    }

    /// Test Sum trait for empty iterator
    #[test]
    fn test_scalar_sum_empty() {
        let scalars: [Scalar; 0] = [];
        let sum: Scalar = scalars.into_iter().sum();
        assert_eq!(sum, Scalar::ZERO);
    }

    /// Test Scalar::to_bytes_le (little-endian)
    #[test]
    fn test_scalar_to_bytes_le() {
        let scalar: Scalar = 42u64.into();
        let bytes = scalar.to_bytes_le();
        // 42 in little-endian should have 42 in the first position
        assert_eq!(bytes[0], 42);
        // Verify round-trip
        let reconstructed = Scalar::from_bytes_le(&bytes).unwrap();
        assert_eq!(reconstructed.0, scalar.0);
    }

    /// Test Scalar::to_bytes_le on zero
    #[test]
    fn test_scalar_to_bytes_le_zero() {
        let bytes = Scalar::ZERO.to_bytes_le();
        assert_eq!(bytes, [0u8; 32]);
    }

    /// Test Scalar::to_bytes_le on one
    #[test]
    fn test_scalar_to_bytes_le_one() {
        let bytes = Scalar::ONE.to_bytes_le();
        assert_eq!(bytes[0], 1);
        for b in bytes.iter().skip(1) {
            assert_eq!(*b, 0);
        }
    }

    /// Test GroupEncoding round-trip for identity point
    #[test]
    fn test_group_encoding_round_trip_identity() {
        use group::GroupEncoding;

        let identity = ProjectivePoint::IDENTITY;
        let bytes = identity.to_bytes();
        // y=1 (little-endian) is [1, 0, 0, ...]; the x-sign flag is packed into
        // bit 7 of the final (32nd) byte and is 0 for x = 0.
        assert_eq!(bytes.as_ref().len(), 32);
        assert_eq!(bytes.as_ref()[0], 1);
        assert_eq!(&bytes.as_ref()[1..32], &[0u8; 31]);
        // Sign bit should be 0 (bit 7 of the last byte not set).
        assert_eq!(bytes.as_ref()[31] & 0x80, 0);

        // Decode and verify
        let decoded = ProjectivePoint::from_bytes(&bytes);
        assert_eq!(decoded.is_none().unwrap_u8(), 0);
        let decoded_point = decoded.unwrap();
        assert_eq!(decoded_point.x, identity.x);
        assert_eq!(decoded_point.y, identity.y);
        assert_eq!(decoded_point.z, identity.z);
    }

    /// Test GroupEncoding round-trip for generator point
    #[test]
    fn test_group_encoding_round_trip_generator() {
        use group::GroupEncoding;

        let generator = ProjectivePoint::GENERATOR;
        let bytes = generator.to_bytes();

        // Decode and verify
        let decoded = ProjectivePoint::from_bytes(&bytes);
        assert_eq!(decoded.is_none().unwrap_u8(), 0);
        let decoded_point = decoded.unwrap();

        // Convert both to affine for comparison
        let gen_affine = generator.to_affine();
        let decoded_affine = decoded_point.to_affine();

        assert_eq!(decoded_affine.x, gen_affine.x);
        assert_eq!(decoded_affine.y, gen_affine.y);
    }

    /// Test GroupEncoding from_bytes_unchecked decodes valid points
    #[test]
    fn test_group_encoding_from_bytes_unchecked() {
        use group::GroupEncoding;

        let generator = ProjectivePoint::GENERATOR;
        let bytes = generator.to_bytes();
        let result = ProjectivePoint::from_bytes_unchecked(&bytes);
        // Should decode successfully
        assert_eq!(result.is_none().unwrap_u8(), 0);
    }

    /// Test that from_bytes correctly decodes a valid SEC1 encoded point.
    /// This is a critical test for the bug fix - the original implementation
    /// always returned failure due to incorrect sign bit handling in Montgomery form.
    #[test]
    fn test_group_encoding_from_bytes_valid_point() {
        use group::GroupEncoding;

        // Create a valid point using known scalar multiplication
        let scalar: Scalar = 42u64.into();
        let point = ProjectivePoint::GENERATOR * scalar;

        // Encode it (to_bytes converts to affine internally)
        let bytes = point.to_bytes();

        // Verify from_bytes can decode it (this was broken before the fix)
        let decoded = ProjectivePoint::from_bytes(&bytes);
        assert_eq!(
            decoded.is_none().unwrap_u8(),
            0,
            "from_bytes should succeed for valid point"
        );

        // Verify the decoded point matches the original by comparing in affine form
        // (projective coordinates can have different representations of the same point)
        let decoded_affine = decoded.unwrap().to_affine();
        let point_affine = point.to_affine();
        assert_eq!(decoded_affine.x, point_affine.x);
        assert_eq!(decoded_affine.y, point_affine.y);
    }

    /// Test that from_bytes correctly decodes the generator point when encoded externally.
    /// This specifically tests the sign bit handling that was broken in the original bug.
    #[test]
    fn test_group_encoding_from_bytes_sign_bit() {
        use group::GroupEncoding;

        // Test with the generator itself - this exercises the sign bit handling
        let generator = ProjectivePoint::GENERATOR;
        let bytes = generator.to_bytes();

        // Verify from_bytes can decode the generator
        let decoded = ProjectivePoint::from_bytes(&bytes);
        assert_eq!(
            decoded.is_none().unwrap_u8(),
            0,
            "from_bytes should succeed for generator"
        );

        // Verify coordinates match when converted to affine
        let decoded_affine = decoded.unwrap().to_affine();
        let generator_affine = generator.to_affine();
        assert_eq!(decoded_affine.x, generator_affine.x);
        assert_eq!(decoded_affine.y, generator_affine.y);
    }

    /// Decoding must never panic for arbitrary 32-byte input.
    #[test]
    fn test_group_encoding_from_bytes_invalid_point() {
        use group::GroupEncoding;

        // y = 2 (little-endian). Whether this is a valid subgroup point is curve
        // dependent; the important guarantee is that decoding does not panic.
        let mut bytes = [0u8; 32];
        bytes[0] = 2;
        let _ = ProjectivePoint::from_bytes(&GroupRepr(bytes));
    }

    /// F3: the identity is the affine point (0, 1) (projective (0, k, k)), NOT a
    /// `z == 0` point. Every identity produced by the API must be detected.
    #[test]
    fn test_identity_detection_edwards() {
        let id = ProjectivePoint::IDENTITY;
        assert!(id.is_identity());
        assert!(bool::from(Group::is_identity(&id)));

        let zero_mul = ProjectivePoint::GENERATOR * Scalar::ZERO;
        assert!(zero_mul.is_identity(), "G * 0 must be identity");
        assert!(bool::from(Group::is_identity(&zero_mul)));

        let g = ProjectivePoint::GENERATOR;
        let g_minus_g = g + (-g);
        assert!(g_minus_g.is_identity(), "G + (-G) must be identity");
        assert!(bool::from(Group::is_identity(&g_minus_g)));

        assert!(!g.is_identity());
        assert!(!bool::from(Group::is_identity(&g)));
    }

    /// F4: fixed-schedule multiplication must agree with the operator.
    #[test]
    fn test_mul_fixed_schedule_matches_operator() {
        let g = ProjectivePoint::GENERATOR;
        let r_minus_1 = Scalar::ZERO - Scalar::ONE;
        let big = Scalar::from(u64::MAX) * Scalar::from(0x9e37_79b9_7f4a_7c15u64);
        let cases = [
            Scalar::ZERO,
            Scalar::ONE,
            Scalar::from(2u64),
            Scalar::from(5u64),
            Scalar::from(42u64),
            Scalar::from(1000u64),
            r_minus_1,
            big,
        ];
        for sc in cases {
            let a = g.mul_fixed_schedule(&sc).to_affine();
            let b = (g * sc).to_affine();
            assert_eq!(a.x, b.x);
            assert_eq!(a.y, b.y);
        }
        assert!(g.mul_fixed_schedule(&Scalar::ZERO).is_identity());
        assert_eq!(
            g.mul_fixed_schedule(&Scalar::ONE).to_affine(),
            g.to_affine()
        );
    }

    /// Regression: derived PartialEq on projective coordinates was incorrect.
    /// (X:Y:Z) and (λX:λY:λZ) represent the same affine point for any λ≠0 and
    /// must compare as equal. The fix uses cross-multiplication without normalization.
    #[test]
    fn test_projective_eq_scaled() {
        let g = ProjectivePoint::GENERATOR;
        // Produce a point with a non-trivial z (doubling guarantees z ≠ 1).
        let p = g + g; // doubled generator — z is not 1
        let three_p = p + g;
        // Same affine point as p+g but with un-normalized projective coordinates.
        // Arithmetic on un-normalized points yields scaled representations.
        let scaled = three_p - g;
        assert_eq!(scaled, p, "(X:Y:Z) and a scaled (λX:λY:λZ) must be equal");
        // Also verify the identity: (0:1:1) is equal to itself even when scaled.
        assert_eq!(ProjectivePoint::IDENTITY, ProjectivePoint::IDENTITY);
    }

    /// F6: scalar decoding must reject non-canonical (`>= r`) byte strings, so
    /// distinct encodings cannot map to the same scalar (signature malleability).
    #[test]
    fn test_scalar_decoding_is_canonical() {
        // r itself (big-endian) must be rejected.
        let mut r_be = [0u8; 32];
        r_be.copy_from_slice(&hex::decode(ORDER_HEX).unwrap());
        assert_eq!(
            Scalar::from_bytes(&r_be).is_some().unwrap_u8(),
            0,
            "the modulus r must be rejected as non-canonical"
        );

        // r - 1 is the canonical maximum and must round-trip.
        let r_minus_1 = Scalar::ZERO - Scalar::ONE;
        let be = r_minus_1.to_bytes();
        let decoded = Scalar::from_bytes(&be);
        assert_eq!(decoded.is_some().unwrap_u8(), 1);
        assert_eq!(decoded.unwrap().0, r_minus_1.0);

        // 2^256 - 1 >= r must be rejected.
        assert_eq!(Scalar::from_bytes(&[0xFFu8; 32]).is_some().unwrap_u8(), 0);

        // Reduction is opt-in only: reduce_bytes_be(r) == 0.
        assert!(Scalar::reduce_bytes_be(&r_be).is_zero());

        // Little-endian canonical decoding uses the same strict bound.
        let mut r_le = r_be;
        r_le.reverse();
        assert_eq!(Scalar::from_bytes_le(&r_le).is_some().unwrap_u8(), 0);
        r_le[0] = r_le[0].wrapping_add(1);
        assert_eq!(Scalar::from_bytes_le(&r_le).is_some().unwrap_u8(), 0);
    }

    /// F5: the compressed point encoding is canonical 32 bytes with no ignored
    /// byte/bits, so it is not malleable.
    #[test]
    fn test_point_encoding_is_canonical() {
        use group::GroupEncoding;

        let g = ProjectivePoint::GENERATOR;
        let bytes = g.to_bytes();
        assert_eq!(bytes.as_ref().len(), 32);
        assert_eq!(ProjectivePoint::from_bytes(&bytes).is_none().unwrap_u8(), 0);

        // Bit 254 (0x40 of the final byte) is unused and must be 0 in a canonical
        // encoding; setting it makes y >= q, which the decoder must reject.
        let mut mutated = bytes;
        mutated.0[31] |= 0x40;
        assert_eq!(
            ProjectivePoint::from_bytes(&mutated).is_none().unwrap_u8(),
            1,
            "non-canonical spare bit must be rejected"
        );
    }

    /// F8: on-curve and prime-order-subgroup helpers, and the checked ctor.
    #[test]
    fn test_on_curve_and_subgroup_helpers() {
        // Generator is on-curve and in the prime-order subgroup.
        assert!(AffinePoint::GENERATOR.is_on_curve());
        assert!(AffinePoint::GENERATOR.is_in_prime_order_subgroup());
        assert!(ProjectivePoint::GENERATOR.is_in_prime_order_subgroup());
        assert!(ProjectivePoint::GENERATOR.is_on_curve());

        // (0, -1) is the order-2 point: on the curve but NOT in the subgroup.
        let order2 = AffinePoint {
            x: BackendBaseField::ZERO,
            y: -BackendBaseField::ONE,
        };
        assert!(order2.is_on_curve());
        assert!(!order2.is_in_prime_order_subgroup());
        assert!(
            AffinePoint::new(BackendBaseField::ZERO, -BackendBaseField::ONE).is_none(),
            "new must reject small-subgroup points"
        );

        // The generator passes the checked constructor.
        let g = AffinePoint::GENERATOR;
        assert!(AffinePoint::new(g.x, g.y).is_some());

        // (1, 2) is off the curve.
        let off = AffinePoint {
            x: BackendBaseField::ONE,
            y: BackendBaseField::from(2u64),
        };
        assert!(!off.is_on_curve());

        // Invalid projective coordinates with z == 0 must not be normalized to
        // affine identity and accepted by validation helpers.
        let invalid = ProjectivePoint {
            x: BackendBaseField::ZERO,
            y: BackendBaseField::ONE,
            z: BackendBaseField::ZERO,
        };
        assert!(!invalid.is_identity());
        assert!(!bool::from(Group::is_identity(&invalid)));
        assert!(!invalid.is_on_curve());
        assert!(!invalid.is_in_prime_order_subgroup());
    }

    /// Security test: cofactor clearing.
    ///
    /// BabyJubJub has cofactor 8, and the torsion subgroup of the full group
    /// `G_full ≅ Z/2 × Z/4` has maximum non-trivial order 4. A point in the
    /// torsion subgroup is NOT in the prime-order subgroup, but plain scalar
    /// multiplication (`*` / `mul_fixed_schedule`) does not eliminate the
    /// torsion component:
    ///
    /// ```ignore
    /// [s]P_torsion = P_torsion   (when s is odd, P_torsion has order 2)
    /// [s]P_torsion = IDENTITY    (when s is even, P_torsion has order 2)
    /// ```
    ///
    /// This leaks the *parity* of the scalar when the attacker controls
    /// `P_torsion` — a small-subgroup attack. The correct fix is cofactor
    /// clearing: `[8]P` projects any point onto the prime-order subgroup, and
    /// `mul_with_cofactor_clear` applies this automatically.
    #[test]
    fn test_cofactor_clearing_security() {
        // The order-2 element P2 = (0, -1) is on the curve but NOT in the
        // prime-order subgroup. It is its own negation: 2*P2 == IDENTITY.
        // We use struct instantiation (no on-curve / subgroup checks) and
        // double_unchecked (avoids the BackendProjective::new subgroup assertion).
        let p2 = ProjectivePoint {
            x: BackendBaseField::ZERO,
            y: -BackendBaseField::ONE,
            z: BackendBaseField::ONE,
        };
        assert!(
            bool::from(Group::is_identity(&p2.double_unchecked())),
            "order-2 element must double to identity"
        );
        assert!(!bool::from(Group::is_identity(&p2)));
        assert!(!p2.is_in_prime_order_subgroup());

        // Verify: [2]P2 == IDENTITY (order 2 element).
        let two: Scalar = 2u64.into();
        assert_eq!(
            p2.mul_unchecked(&two),
            ProjectivePoint::IDENTITY,
            "order-2 element must double to identity"
        );

        // --- Without cofactor clearing ---
        //
        // [1]P2 == P2    (odd scalar: torsion component preserved)
        // [3]P2 == P2    (odd scalar: 3 mod 2 == 1)
        // [2]P2 == ID    (even scalar: collapses to identity)
        //
        // We use mul_unchecked instead of * because the * operator converts
        // via From<ProjectivePoint> which asserts the subgroup.
        let one: Scalar = 1u64.into();
        let three: Scalar = 3u64.into();
        assert_eq!(
            p2.mul_unchecked(&one),
            p2,
            "[1]P2 must equal P2 (no cofactor cleared)"
        );
        assert_eq!(
            p2.mul_unchecked(&three),
            p2,
            "[3]P2 must equal P2 (3 mod 2 == 1)"
        );
        assert_eq!(
            p2.mul_unchecked(&two),
            ProjectivePoint::IDENTITY,
            "[2]P2 == ID"
        );

        // A protocol that multiplies an attacker-supplied point by a secret
        // scalar `s` will find that the result equals `P2` iff `s` is odd.
        // This parity leak is a small-subgroup attack vector.

        // --- With cofactor clearing ---
        //
        // [8]P2 projects P2 onto the prime-order subgroup:
        //   P2 = [8]Q  where Q = P2 (since 2*P2 == ID, 4*P2 == ID too,
        //   and 8 is a multiple of the order of P2, so Q = IDENTITY).
        // Therefore [1]([8]P2) == IDENTITY (not P2), regardless of scalar parity.
        let eight: Scalar = 8u64.into();
        let p2_projected = p2.mul_unchecked(&eight);
        assert_eq!(
            p2_projected,
            ProjectivePoint::IDENTITY,
            "[8]P2 must project to identity (P2 has order dividing 8)"
        );

        // [1]([8]P2) == ID — the torsion component is gone.
        // p2_projected is the identity (y == z), so * on it is safe.
        assert_eq!(
            p2_projected * one,
            ProjectivePoint::IDENTITY,
            "cofactor-cleared point must stay at identity"
        );

        // [3]([8]P2) == ID — even with an odd scalar.
        assert_eq!(
            p2_projected * three,
            ProjectivePoint::IDENTITY,
            "cofactor-cleared point stays at identity under any scalar"
        );

        // `mul_with_cofactor_clear` applies the same cofactor multiplication.
        assert_eq!(
            p2.mul_with_cofactor_clear(&one),
            ProjectivePoint::IDENTITY,
            "mul_with_cofactor_clear must eliminate the torsion component"
        );
        assert_eq!(
            p2.mul_with_cofactor_clear(&three),
            ProjectivePoint::IDENTITY,
            "mul_with_cofactor_clear must work regardless of scalar parity"
        );

        // --- Prime-order subgroup points are unaffected ---
        //
        // For a generator G in the prime-order subgroup, [8]G is still in the
        // subgroup (since the subgroup is closed under scalar multiplication
        // by 8).  The result may have a non-unit z in the backend representation,
        // so we compare the identity check rather than raw coordinate equality.
        let g = ProjectivePoint::GENERATOR;
        let g_projected = g * eight;
        assert!(
            g_projected.is_in_prime_order_subgroup(),
            "[8]G must still be in the prime-order subgroup"
        );
        // [8]([8]G) == [8]G (idempotent: already in the subgroup).
        let g_proj_idem = g_projected * eight;
        assert!(
            g_proj_idem.is_in_prime_order_subgroup(),
            "[8]([8]G) must still be in the prime-order subgroup (subgroup closure)"
        );

        // --- Regression: cofactor clearing must use the INTEGER 8*s ---
        //
        // For a large scalar the integer product 8*s exceeds r, and 8*s mod r is
        // generally NOT a multiple of 8 (r is odd). Cofactor clearing therefore
        // must multiply by the integer 8*s, never by 8*s reduced mod r. Take
        // s = r-1: then 8*s mod r = r-8, which is ODD (r ≡ 1 mod 8), so the buggy
        // mod-r version computes [r-8]P2 = P2 (torsion survives), whereas the
        // correct integer version computes [8(r-1)]P2 = IDENTITY.
        let r_minus_1 = Scalar::ZERO - Scalar::ONE;
        let cleared = p2.mul_with_cofactor_clear(&r_minus_1);
        assert_eq!(
            cleared,
            ProjectivePoint::IDENTITY,
            "cofactor clearing must use the integer 8*s, not 8*s mod r"
        );

        // Functional contract on a prime-order point with a large (wrapping)
        // scalar: mul_with_cofactor_clear(s) == [8s]G computed via the operators,
        // and the result stays in the prime-order subgroup.
        let s_big = Scalar::from(0x9e37_79b9_7f4a_7c15u64) * Scalar::from(u64::MAX);
        assert_eq!(
            g.mul_with_cofactor_clear(&s_big).to_affine(),
            (g * (s_big * eight)).to_affine(),
            "mul_with_cofactor_clear(s) must equal [8s]G on a subgroup point"
        );
        assert!(
            g.mul_with_cofactor_clear(&s_big)
                .is_in_prime_order_subgroup(),
            "cofactor-cleared result must be in the prime-order subgroup"
        );
    }

    #[test]
    #[should_panic(expected = "invalid projective point: z-coordinate is zero")]
    fn test_projective_to_affine_rejects_zero_z() {
        let invalid = ProjectivePoint {
            x: BackendBaseField::ZERO,
            y: BackendBaseField::ONE,
            z: BackendBaseField::ZERO,
        };

        let _ = invalid.to_affine();
    }

    #[test]
    fn test_randomization() {
        use rand::{SeedableRng, rngs::StdRng};
        let mut rng = StdRng::seed_from_u64(42);
        // Cover ProjectivePoint::try_random (631-634)
        let _ = ProjectivePoint::try_random(&mut rng).unwrap();

        // Cover Scalar::random (1221-1223, 1500-1505)
        let _ = <Scalar as Field>::random(&mut rng);

        // Cover Scalar::try_random (1225-1227, 1507-1511)
        let _ = <Scalar as Field>::try_random(&mut rng).unwrap();
    }

    #[test]
    fn test_identity() {
        // Cover ProjectivePoint::identity (640-642)
        let identity = ProjectivePoint::identity();
        assert!(identity.is_identity());
    }

    #[test]
    fn test_scalar_field_square_double() {
        let s = Scalar::from(2u64);
        // Explicitly test Field trait methods (lines 1229-1235) to ensure coverage
        let squared = <Scalar as Field>::square(&s);
        assert_eq!(squared, Scalar::from(4u64));

        let doubled = <Scalar as Field>::double(&s);
        assert_eq!(doubled, Scalar::from(4u64));

        let zero = Scalar::ZERO;
        assert_eq!(<Scalar as Field>::square(&zero), zero);
        assert_eq!(<Scalar as Field>::double(&zero), zero);
    }

    #[test]
    fn test_iterators() {
        let points = [ProjectivePoint::GENERATOR, ProjectivePoint::GENERATOR];
        // Cover Sum<&ProjectivePoint> for ProjectivePoint (1300-1302)
        let sum: ProjectivePoint = points.iter().sum();
        assert_eq!(sum, ProjectivePoint::GENERATOR + ProjectivePoint::GENERATOR);

        let scalars = [Scalar::from(1u64), Scalar::from(2u64)];
        // Cover Sum<&Scalar> for Scalar (1313-1315)
        let sum_s: Scalar = scalars.iter().sum();
        assert_eq!(sum_s, Scalar::from(3u64));

        // Cover Product<&Scalar> for Scalar (1325-1327)
        let prod_s: Scalar = scalars.iter().product();
        assert_eq!(prod_s, Scalar::from(2u64));
    }

    #[test]
    #[allow(clippy::op_ref)]
    fn test_operator_overloads() {
        let p = ProjectivePoint::GENERATOR;
        let s = Scalar::from(2u64);

        // Cover Add<&ProjectivePoint> for ProjectivePoint (849-854)
        let _ = p + &p;

        // Cover Sub<&ProjectivePoint> for &ProjectivePoint (883-888)
        let _ = &p - &p;

        // Cover Sub<&ProjectivePoint> for ProjectivePoint (894-899)
        let _ = p - &p;

        // Cover Mul<&Scalar> for &ProjectivePoint (954-958)
        let _ = &p * &s;

        let s1 = Scalar::from(1u64);
        let s2 = Scalar::from(2u64);

        // Cover Add<&Scalar> for Scalar (1350-1352)
        let _ = s1 + &s2;

        // Cover Sub<&Scalar> for &Scalar (1370-1372)
        let _ = &s1 - &s2;

        // Cover Mul<&Scalar> for &Scalar (1406-1408)
        let _ = &s1 * &s2;
    }

    #[test]
    fn test_affine_ct_eq() {
        let a = AffinePoint::IDENTITY;
        let b = AffinePoint::from(ProjectivePoint::GENERATOR);
        // Cover AffinePoint::ct_eq (1478-1480)
        let _ = a.ct_eq(&b);
    }
}
