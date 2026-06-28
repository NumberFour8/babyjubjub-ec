# babyjubjub-ec

[![Crates.io](https://img.shields.io/crates/v/babyjubjub-ec)](https://crates.io/crates/babyjubjub-ec)
[![Docs](https://docs.rs/babyjubjub-ec/badge.svg)](https://docs.rs/babyjubjub-ec)

BabyJubJub elliptic curve implementation wrapped in `elliptic-curve` traits.

This crate provides a wrapper around the `taceo-ark-babyjubjub` crate that implements
the BabyJubJub curve in a way compatible with the `elliptic-curve` crate traits.

## What is BabyJubJub?

BabyJubJub is a twisted Edwards curve that is birationally equivalent to the Edwards curve
ed25519. It is designed to be efficient for arithmetic circuits and is commonly used in
zero-knowledge proofs like zk-SNARKs and ZK-Rollups.

- **Curve type**: Twisted Edwards curve (cofactor 8)
- **Base field `Fq`**: a **prime** field, where `q` is a **254-bit** prime (the BN254 scalar field). Point coordinates live in `Fq`.
- **Scalar field `Fr`** (`Scalar`): a **prime** field of **251-bit** prime order `r` (the prime-order subgroup):
  `r = 2736030358979909402780800718157159386076813972158567259200215660948447373041`

## Features

- `std` (default): enables standard-library support for this wrapper and the
  arkworks crates it uses directly.
- `zeroize` (**off** by default): retained for backwards compatibility, but now a
  no-op. `zeroize::DefaultIsZeroes` (and hence an explicit `.zeroize()`) is
  implemented for `Scalar`, `AffinePoint`, and `ProjectivePoint`
  **unconditionally** via `elliptic-curve`'s re-exported `zeroize`, which
  `CurveArithmetic` requires — so it no longer depends on this feature.
- Disabling default features builds `no_std`. Note that the arkworks backend still
  requires a global allocator (`alloc`), so this targets `no_std + alloc`
  environments rather than bare-metal `no_std` without an allocator.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
babyjubjub-ec = "0.3"
```

## Usage

### Scalar Multiplication

```rust
use babyjubjub_ec::{ProjectivePoint, Scalar};

// Create a scalar from a u64
let scalar: Scalar = 42u64.into();

// Multiply the generator by the scalar
let result = ProjectivePoint::GENERATOR * scalar;

// Convert to affine coordinates for reading
let affine = result.to_affine();
println!("Result x: {:?}", affine.x);
println!("Result y: {:?}", affine.y);
```

### Point Addition

```rust
use babyjubjub_ec::{ProjectivePoint, Scalar};

// Get the generator point
let g = ProjectivePoint::GENERATOR;

// Compute 2 * G using addition
let two_g = g + g;

// Compute 3 * G using addition
let three_g = two_g + g;

// Verify: 3 * G = G + G + G
let expected = ProjectivePoint::GENERATOR * Scalar::from(3u64);
assert_eq!(three_g.to_affine().x, expected.to_affine().x);
assert_eq!(three_g.to_affine().y, expected.to_affine().y);
```

### Working with Scalars

```rust
use babyjubjub_ec::{Scalar, ProjectivePoint};

// Create scalars
let a: Scalar = 10u64.into();
let b: Scalar = 20u64.into();

// Scalar arithmetic
let sum = a + b;        // 10 + 20 = 30
let product = a * b;    // 10 * 20 = 200

// Scalar multiplication with generator
let result = ProjectivePoint::GENERATOR * sum;
```

### Serialization

```rust
use babyjubjub_ec::{ProjectivePoint, GroupRepr};

let point = ProjectivePoint::GENERATOR;

// Serialize to the canonical 32-byte compressed encoding
// (little-endian y with the x-sign flag packed into the top bit of the last byte).
let bytes: GroupRepr = point.to_bytes();

// Deserialize from bytes. `from_bytes` returns a `CtOption` and validates that
// the point is on-curve AND in the prime-order subgroup; non-canonical or
// off-curve/small-subgroup encodings are rejected.
let decoded = ProjectivePoint::from_bytes(&bytes);
```

## API Overview

### Types

| Type | Description |
|------|-------------|
| `BabyJubJub` | Curve type implementing `Curve`, `PrimeCurve`, and `CurveArithmetic` |
| `AffinePoint` | Affine point representation (x, y coordinates) |
| `ProjectivePoint` | Projective point representation (x, y, z coordinates) |
| `Scalar` | Scalar field element |
| `GroupRepr` | Group element representation for serialization |
| `FieldElement` | Base-field intermediate used by hash-to-curve (`GroupDigest::FieldElement`) |

### Traits Implemented

- `elliptic_curve::Curve` and `elliptic_curve::PrimeCurve` for `BabyJubJub`
- `elliptic_curve::CurveArithmetic` for `BabyJubJub`, with `AffinePoint`,
  `ProjectivePoint`, and `Scalar` as its associated types. This requires (and the
  crate provides) the full set of helper-trait impls, including:
  - `group::Curve` and `BatchNormalize` for `ProjectivePoint`
  - `elliptic_curve::point::AffineCoordinates` for `AffinePoint` (the identity /
    generator / `to_curve` / `from_coordinates` constructors and the `y` /
    `x_is_odd` getters are inherent methods, since `elliptic-curve` 0.13.8 has no
    `CurveAffine`)
  - `elliptic_curve::ops::{Invert, Reduce, MulByGenerator, LinearCombination}`,
    `elliptic_curve::scalar::{FromUintUnchecked, IsHigh}`, `core::ops::ShrAssign`,
    and `bigint::modular::Retrieve` for `Scalar`
  - conversions to/from `U256`, `ScalarPrimitive`, `FieldBytes`, `NonZeroScalar`,
    and `NonIdentity`, plus `Scalar * Point` multiplication in both operand orders
- `group::Group` for `ProjectivePoint`
- `group::prime::PrimeGroup` and `group::cofactor::CofactorGroup` for
  `ProjectivePoint` (with `Subgroup = ProjectivePoint`). `clear_cofactor`
  multiplies by the cofactor 8 to map any point — including torsion points — into
  the prime-order subgroup, while `is_torsion_free` / `into_subgroup` report
  subgroup membership.
- `group::ff::Field` for `Scalar`
- `group::ff::PrimeField` for `Scalar`. **Note:** `PrimeField::Repr` is now
  `elliptic_curve::FieldBytes<BabyJubJub>` (i.e. `Array<u8, U32>`) rather than
  `[u8; 32]`, as `CurveArithmetic` mandates. If you have a `[u8; 32]`, pass
  `bytes.into()` to `from_repr` / `from_repr_vartime`.
- `group::GroupEncoding` for `ProjectivePoint` and `AffinePoint`
- `elliptic_curve::hash2curve::GroupDigest` for `BabyJubJub` (RFC 9380
  hash-to-curve). `hash_from_bytes` / `encode_from_bytes` map arbitrary byte
  strings to prime-order-subgroup points via the twisted-Edwards **Elligator2**
  map (`FieldElement` as `GroupDigest::FieldElement`) followed by cofactor
  clearing; `hash_to_scalar` is available since `Scalar` implements
  `hash2curve::FromOkm`. The expander/hash is caller-chosen, e.g.
  `BabyJubJub::hash_from_bytes::<ExpandMsgXmd<sha2::Sha256>>(&[msg], &[dst])`.
- `subtle::ConditionallySelectable` and `subtle::ConstantTimeEq` for
  `ProjectivePoint`, `AffinePoint`, `Scalar`
- `zeroize::DefaultIsZeroes` for all point and scalar types (**unconditional**)

### `elliptic-curve` 0.13 line

This is the **`legacy-ec-13`** release line (crate version `0.2.x`), built against
the **`elliptic-curve` 0.13.8** ecosystem (`ff` 0.13 / `group` 0.13 /
`rand_core` 0.6). It provides the same `CurveArithmetic` functionality as the
`elliptic-curve` 0.14 line (crate `0.3.x`), re-expressed with 0.13.8 conventions.
APIs that exist only in 0.14 are therefore **not** available here; their
behaviour stays reachable through the 0.13.8-native equivalents:

- `ctutils::{CtEq, CtSelect}` → use `subtle::{ConstantTimeEq, ConditionallySelectable}`.
- `Generate` / `try_random` → use `group::Group::random` / `ff::Field::random`.
- `MulVartime` / `MulByGeneratorVartime` → use the constant-time `*` operator and
  `MulByGenerator`.
- `CurveAffine` → identity / generator / `to_curve` are inherent on `AffinePoint`.
- the `ScalarValue` type → named `ScalarPrimitive` in 0.13.8.

## Examples

See the [tests](src/lib.rs) for more examples including:

- Scalar arithmetic (add, sub, mul, invert, square)
- Point operations (add, sub, neg, double)
- Conversion between affine and projective coordinates
- Random point generation
- Fixed-schedule scalar multiplication and constant-time selection/equality operations

## Testing

```bash
cargo test
cargo test --all-features
cargo check --no-default-features
```

## Benchmarking

```bash
cargo bench
```

## Security

This crate is a thin wrapper over the arkworks backend. Please note:

- **Scalar multiplication is *almost* constant-time.** Both the `*` operator and
  `ProjectivePoint::mul_fixed_schedule` run a constant-time scalar-multiplication
  *algorithm* (no scalar-dependent loop length, branching, or memory-access
  pattern), but through two distinct code paths. The `*` operator delegates to the
  backend's `mul_projective`, which `taceo-ark-babyjubjub` implements as a
  **Montgomery ladder** over a fixed number of scalar bits (no leading-zero
  skipping) with branch-free, bit-masked register swaps.
  `ProjectivePoint::mul_fixed_schedule`, by contrast, does **not** use the
  backend's `mul_projective`: it runs an in-crate **double-and-add-always** loop
  built only on the curve's complete point addition and doubling, so its
  algorithm-level constant-time property does not depend on the backend's
  scalar-multiplication routine. Neither is **end-to-end** constant-time, however:
  the underlying `ark-ff` field arithmetic uses a data-dependent conditional
  reduction (`Fp::subtract_modulus`), leaving a small residual timing signal. For
  end-to-end constant time, use a backend with bit-masked field reduction
  throughout.
- **Variable-time scalar-field operations.** `Scalar::invert` and `Scalar`'s
  `sqrt`/`sqrt_ratio` delegate to the backend and are **not** constant-time
  (input-dependent control flow).
- **Validation.** `ProjectivePoint::from_bytes` validates on-curve and
  prime-order-subgroup membership. The `from_bytes_unchecked` decoder does **not**;
  for untrusted coordinates use `AffinePoint::new` or the `is_on_curve` /
  `is_in_prime_order_subgroup` helpers.
- **Canonical encodings.** Scalar decoding (`from_bytes`, `from_repr`) rejects
  non-canonical values `>= r`, and the point encoding is a canonical 32 bytes,
  preventing scalar/point malleability. Use `Scalar::reduce_bytes_be` when
  modular reduction is explicitly desired.
- **Zeroization.** `Scalar` is `Copy`, so it cannot auto-zeroize on drop; wipe
  secret storage yourself. The point and scalar types implement `Zeroize` via
  `DefaultIsZeroes` **unconditionally** (required by `CurveArithmetic`).

## License

GPL-3.0
