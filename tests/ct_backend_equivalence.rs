//! Cross-backend equivalence (constant-time `fiat` backend vs. reference).
//!
//! These tests are compiled only with `--features fiat`. They drive the crate's
//! public API — which, under `fiat`, runs entirely over the constant-time
//! fiat-crypto field/curve backend — and compare its results, via the canonical
//! 32-byte point encoding, against an **independent** computation done directly
//! with the default `taceo-ark-babyjubjub` reference implementation.
//!
//! This is an end-to-end oracle for the M1/M2 mitigation: scalar multiplication
//! (`*`, `mul_fixed_schedule`, `mul_with_cofactor_clear`), point decoding (whose
//! decompression performs an `Fq` square root) and the generator must all match
//! the reference bit-for-bit, proving the constant-time backend is behaviorally
//! identical to the fast one.
#![cfg(feature = "fiat")]

use ark_ec::{CurveGroup, PrimeGroup};
use ark_serialize::{CanonicalSerialize, Compress};
use babyjubjub_ec::{GroupRepr, ProjectivePoint, Scalar};
use group::GroupEncoding;
use taceo_ark_babyjubjub::{EdwardsProjective as TaceoProj, Fr as TaceoFr};

/// Canonical 32-byte compressed encoding of a taceo projective point (the same
/// format the crate emits from `to_bytes`).
fn taceo_encode(p: TaceoProj) -> [u8; 32] {
    let affine = p.into_affine();
    let mut out = [0u8; 32];
    affine
        .serialize_with_mode(&mut out[..], Compress::Yes)
        .expect("compressed serialization is 32 bytes");
    out
}

/// Deterministic SplitMix64, to derive reproducible scalars without an rng dep.
fn next_u64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Build the *same* scalar value in both representations from three `u64`s,
/// using only field arithmetic (`a*b + c`). This avoids any byte-endianness
/// assumptions while still producing full-width (≈251-bit) scalars.
fn paired_scalar(a: u64, b: u64, c: u64) -> (Scalar, TaceoFr) {
    let crate_s = Scalar::from(a) * Scalar::from(b) + Scalar::from(c);
    let taceo_s = TaceoFr::from(a) * TaceoFr::from(b) + TaceoFr::from(c);
    (crate_s, taceo_s)
}

#[test]
fn generator_encoding_matches_reference() {
    let crate_g = ProjectivePoint::GENERATOR.to_bytes();
    let taceo_g = taceo_encode(TaceoProj::generator());
    assert_eq!(
        crate_g.as_ref(),
        taceo_g.as_slice(),
        "crate generator must encode identically to the taceo reference"
    );
}

#[test]
fn small_scalar_mul_matches_reference() {
    for k in [0u64, 1, 2, 3, 7, 8, 9, 1000, u32::MAX as u64, u64::MAX] {
        let got = (ProjectivePoint::GENERATOR * Scalar::from(k)).to_bytes();
        let expected = taceo_encode(TaceoProj::generator() * TaceoFr::from(k));
        assert_eq!(
            got.as_ref(),
            expected.as_slice(),
            "[{k}]G disagreed with the taceo reference"
        );
    }
}

#[test]
fn scalar_mul_matches_reference() {
    let mut st = 0xbaba_d1ce_5eed_0001u64;
    for _ in 0..48 {
        let (s, ts) = paired_scalar(next_u64(&mut st), next_u64(&mut st), next_u64(&mut st));
        let got = (ProjectivePoint::GENERATOR * s).to_bytes();
        let expected = taceo_encode(TaceoProj::generator() * ts);
        assert_eq!(
            got.as_ref(),
            expected.as_slice(),
            "operator scalar-mul disagreed with the taceo reference"
        );
    }
}

#[test]
fn mul_fixed_schedule_matches_reference() {
    // `mul_fixed_schedule` is the in-crate double-and-add-always routine; it is
    // markedly slower than the ladder on the (unoptimized) fiat backend, so use
    // a smaller sample — equivalence is structural, not statistical.
    let mut st = 0xfeed_5ced_0000_0001u64;
    for _ in 0..8 {
        let (s, ts) = paired_scalar(next_u64(&mut st), next_u64(&mut st), next_u64(&mut st));
        let got = ProjectivePoint::GENERATOR.mul_fixed_schedule(&s).to_bytes();
        let expected = taceo_encode(TaceoProj::generator() * ts);
        assert_eq!(
            got.as_ref(),
            expected.as_slice(),
            "mul_fixed_schedule disagreed with the taceo reference"
        );
    }
}

#[test]
fn mul_with_cofactor_clear_matches_reference() {
    // On a prime-order-subgroup point (the generator), clearing the cofactor is
    // `[8·s]G`; the integer factor 8 is applied without reducing mod r, but G
    // has order r, so the reference can multiply by `Fr::from(8)·s`.
    let mut st = 0xC0FA_C708_0000_0001u64;
    let eight = TaceoFr::from(8u64);
    for _ in 0..16 {
        let (s, ts) = paired_scalar(next_u64(&mut st), next_u64(&mut st), next_u64(&mut st));
        let got = ProjectivePoint::GENERATOR
            .mul_with_cofactor_clear(&s)
            .to_bytes();
        let expected = taceo_encode(TaceoProj::generator() * (eight * ts));
        assert_eq!(
            got.as_ref(),
            expected.as_slice(),
            "mul_with_cofactor_clear disagreed with the taceo reference"
        );
    }
}

#[test]
fn decode_of_reference_encoding_matches() {
    // Decode a taceo-produced encoding with the crate's (fiat) `from_bytes`,
    // which performs an `Fq` square-root decompression, and confirm it yields
    // the same point the crate computes for `[s]G`.
    let mut st = 0xDEC0_DE00_0000_0001u64;
    for _ in 0..16 {
        let (s, ts) = paired_scalar(next_u64(&mut st), next_u64(&mut st), next_u64(&mut st));
        let reference_bytes = taceo_encode(TaceoProj::generator() * ts);

        let decoded = ProjectivePoint::from_bytes(&GroupRepr(reference_bytes));
        assert!(
            bool::from(decoded.is_some()),
            "reference encoding must decode as a valid subgroup point"
        );
        let decoded = decoded.unwrap();

        let direct = ProjectivePoint::GENERATOR * s;
        assert_eq!(
            decoded.to_bytes().as_ref(),
            direct.to_bytes().as_ref(),
            "decoded reference point must equal the directly computed [s]G"
        );
    }
}
