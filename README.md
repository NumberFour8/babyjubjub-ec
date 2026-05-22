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

- **Curve type**: Twisted Edwards curve
- **Prime order**: 2736030358979909402780800718157159386076813972158567259200215660948447373041
- **Field**: Binary extension field Fq where q is a 255-bit prime

## Features

- `std` (default): Enable standard library support
- `no_std`: Available for embedded environments

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

// Serialize to bytes (33 bytes: 32 for y-coordinate + 1 for sign)
let bytes: GroupRepr = point.to_bytes();

// Deserialize from bytes
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
- Constant-time operations

## Testing

```bash
cargo test
cargo test --features std
```

## Benchmarking

```bash
cargo bench
```

## License

GPL-3.0
