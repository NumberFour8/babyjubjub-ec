//! Curve/field backend selection.
//!
//! `lib.rs` builds its public `Scalar` / `AffinePoint` / `ProjectivePoint`
//! types on four backend aliases (`BackendBaseField`, `BackendScalar`,
//! `BackendAffine`, `BackendProjective`) plus the curve generator coordinates
//! (`GENERATOR_X`, `GENERATOR_Y`). This module picks which implementation those
//! names resolve to:
//!
//! - **default** (no `fiat` feature): the fast `taceo-ark-babyjubjub` backend,
//!   whose field arithmetic uses ark-ff's data-dependent conditional reduction
//!   (`subtract_modulus`).
//! - **`fiat` feature**: a constant-time backend whose field arithmetic is the
//!   vendored, formally-verified fiat-crypto code (branch-free masked reduction
//!   plus Bernstein-Yang inversion), with arkworks' twisted-Edwards formulas and
//!   taceo's constant-time Montgomery ladder reused over the new field types.
//!
//! Both backends describe the *same* curve with identical canonical encodings,
//! so the rest of the crate stays generic over these aliases and is unchanged.

#[cfg(feature = "fiat")]
pub mod ct_curve;
#[cfg(feature = "fiat")]
pub mod ct_field;
#[cfg(feature = "fiat")]
pub mod fiat;

// ----- default backend: fast, variable-time field reduction -----
#[cfg(not(feature = "fiat"))]
pub use taceo_ark_babyjubjub::{
    EdwardsAffine as BackendAffine, EdwardsProjective as BackendProjective, Fq as BackendBaseField,
    Fr as BackendScalar, GENERATOR_X, GENERATOR_Y,
};

// ----- `fiat` backend: constant-time fiat-crypto field reduction -----
#[cfg(feature = "fiat")]
pub use ct_curve::{
    CtAffine as BackendAffine, CtProjective as BackendProjective, GENERATOR_X, GENERATOR_Y,
};
#[cfg(feature = "fiat")]
pub use ct_field::{CtFq as BackendBaseField, CtFr as BackendScalar};

// Wrap raw Montgomery-form limbs as a backend field element, without range or
// canonicality checks. This is the backend-agnostic equivalent of ark-ff's
// `Fp::new_unchecked` (which exists only on the `MontBackend` field type): both
// the default and `fiat` field representations store reduced Montgomery limbs
// with `R = 2^256`, so the same hardcoded limb constants are valid for either.
// Used by `lib.rs` to define its compile-time field/scalar constants.
#[cfg(not(feature = "fiat"))]
pub const fn new_base_field_unchecked(limbs: ark_ff::BigInt<4>) -> BackendBaseField {
    BackendBaseField::new_unchecked(limbs)
}
#[cfg(not(feature = "fiat"))]
pub const fn new_scalar_unchecked(limbs: ark_ff::BigInt<4>) -> BackendScalar {
    BackendScalar::new_unchecked(limbs)
}
#[cfg(feature = "fiat")]
pub const fn new_base_field_unchecked(limbs: ark_ff::BigInt<4>) -> BackendBaseField {
    ct_field::ct_fq(limbs)
}
#[cfg(feature = "fiat")]
pub const fn new_scalar_unchecked(limbs: ark_ff::BigInt<4>) -> BackendScalar {
    ct_field::ct_fr(limbs)
}
