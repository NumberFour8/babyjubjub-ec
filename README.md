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

- `std` (default): enables standard-library support (and `std` on the arkworks backend).
- Disabling default features builds `no_std`. Note that the arkworks backend still
  requires a global allocator (`alloc`), so this targets `no_std + alloc`
  environments rather than bare-metal `no_std` without an allocator.

## Installation

Add this to your `Cargo.toml`:

```toml
[dependencies]
babyjubjub-ec = "0.1"
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
| `BabyJubJub` | Curve type implementing `Curve` and `PrimeCurve` |
| `AffinePoint` | Affine point representation (x, y coordinates) |
| `ProjectivePoint` | Projective point representation (x, y, z coordinates) |
| `Scalar` | Scalar field element |
| `GroupRepr` | Group element representation for serialization |

### Traits Implemented

- `group::Group` for `ProjectivePoint`
- `group::ff::Field` for `Scalar`
- `group::ff::PrimeField` for `Scalar`
- `elliptic_curve::GroupEncoding` for `ProjectivePoint`
- `subtle::ConditionallySelectable` for `ProjectivePoint`, `AffinePoint`, `Scalar`
- `zeroize::DefaultIsZeroes` for all point and scalar types

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
cargo test --features std
```

## Benchmarking

```bash
cargo bench
```

## Security

This crate is a thin wrapper over the arkworks backend. Please note:

- **Variable-time arithmetic.** Scalar multiplication via the `*` operator,
  `Scalar::invert`, and `Scalar`'s `sqrt`/`sqrt_ratio` delegate to the backend
  and are **not** constant-time. `ProjectivePoint::mul_fixed_schedule` avoids
  scalar-dependent control flow in this wrapper, but still calls backend group
  operations and must not be treated as an end-to-end constant-time primitive.
- **Validation.** `ProjectivePoint::from_bytes` validates on-curve and
  prime-order-subgroup membership. The raw constructors `AffinePoint::new_unchecked` /
  `ProjectivePoint::new_unchecked` and `from_bytes_unchecked` do **not**; for untrusted
  coordinates use `AffinePoint::new` or the `is_on_curve` /
  `is_in_prime_order_subgroup` helpers.
- **Canonical encodings.** Scalar decoding (`from_bytes`, `from_repr`) rejects
  non-canonical values `>= r`, and the point encoding is a canonical 32 bytes,
  preventing scalar/point malleability. Use `Scalar::reduce_bytes_be` when
  modular reduction is explicitly desired.
- **Zeroization.** `Scalar` is `Copy`, so it cannot auto-zeroize on drop; wipe
  secret storage yourself (the type implements `Zeroize` via `DefaultIsZeroes`).

## License

GPL-3.0
