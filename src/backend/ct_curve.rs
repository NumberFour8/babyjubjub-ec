//! BabyJubJub twisted-Edwards curve configuration over the constant-time field.
//! Compiled only with `--features fiat`.
//!
//! This is a faithful fork of taceo's `EdwardsConfig` (`taceo-ark-babyjubjub`),
//! re-parameterized over the constant-time [`CtFq`]/[`CtFr`] field types so that
//! arkworks' audited, generic twisted-Edwards `add`/`double` formulas execute
//! over branch-free fiat-crypto arithmetic. It keeps taceo's constant-time
//! Montgomery-ladder `mul_projective`/`mul_affine` overrides verbatim.
//!
//! Every curve constant is *reinterpreted* from the corresponding taceo
//! constant (identical `R = 2^256` Montgomery limbs), so the two backends
//! describe exactly the same curve. `MontFp!` cannot be used here because it
//! only supports the `MontBackend` config; the limb-copy via [`ct_fq`]/[`ct_fr`]
//! achieves the same result for the custom constant-time config.

use ark_ec::{
    AdditiveGroup,
    models::CurveConfig,
    twisted_edwards::{Affine, MontCurveConfig, Projective, TECurveConfig},
};
use ark_ff::{BigInt, Zero};
use taceo_ark_babyjubjub::EdwardsConfig as TaceoEdwards;

use super::ct_field::{CtFq, CtFr, ct_fq, ct_fr};

/// Constant-time affine point on BabyJubJub.
pub type CtAffine = Affine<CtEdwardsConfig>;
/// Constant-time projective (extended twisted-Edwards) point on BabyJubJub.
pub type CtProjective = Projective<CtEdwardsConfig>;

/// BabyJubJub twisted-Edwards configuration over the constant-time field.
///
/// Curve equation: `168700·x² + y² = 1 + 168696·x²·y²` over `Fq`.
#[derive(Clone, Default, PartialEq, Eq)]
pub struct CtEdwardsConfig;

impl CurveConfig for CtEdwardsConfig {
    type BaseField = CtFq;
    type ScalarField = CtFr;

    /// COFACTOR = 8 (shared with the default backend).
    const COFACTOR: &'static [u64] = <TaceoEdwards as CurveConfig>::COFACTOR;

    /// COFACTOR^(-1) mod r.
    const COFACTOR_INV: CtFr = ct_fr(<TaceoEdwards as CurveConfig>::COFACTOR_INV.0);
}

impl TECurveConfig for CtEdwardsConfig {
    /// COEFF_A = 168700
    const COEFF_A: CtFq = ct_fq(<TaceoEdwards as TECurveConfig>::COEFF_A.0);

    /// COEFF_D = 168696
    const COEFF_D: CtFq = ct_fq(<TaceoEdwards as TECurveConfig>::COEFF_D.0);

    /// AFFINE_GENERATOR_COEFFS = (GENERATOR_X, GENERATOR_Y)
    const GENERATOR: CtAffine = CtAffine::new_unchecked(GENERATOR_X, GENERATOR_Y);

    type MontCurveConfig = CtEdwardsConfig;

    /// Constant-time scalar multiplication (Montgomery ladder with bit-masked
    /// conditional swaps). Identical to taceo's override: a fixed number of
    /// operations regardless of the scalar bits, so it does not leak the scalar.
    fn mul_projective(base: &Projective<Self>, scalar: &[u64]) -> Projective<Self> {
        let mut r0 = Projective::<Self>::zero();
        let mut r1 = *base;
        let mut prev_bit = false;
        for b in ark_ff::BitIteratorBE::new(scalar) {
            let swap = prev_bit ^ b;
            prev_bit = b;
            conditional_swap(&mut r0, &mut r1, swap);
            r1 += r0;
            r0.double_in_place();
        }
        conditional_select(&mut r0, &r1, prev_bit);
        r0
    }

    /// Also override `mul_affine` to use the constant-time `mul_projective`.
    fn mul_affine(base: &Affine<Self>, scalar: &[u64]) -> Projective<Self> {
        let base = Projective::<Self>::from(*base);
        Self::mul_projective(&base, scalar)
    }
}

impl MontCurveConfig for CtEdwardsConfig {
    /// COEFF_A = 168698
    const COEFF_A: CtFq = ct_fq(<TaceoEdwards as MontCurveConfig>::COEFF_A.0);
    /// COEFF_B = 1
    const COEFF_B: CtFq = ct_fq(<TaceoEdwards as MontCurveConfig>::COEFF_B.0);

    type TECurveConfig = CtEdwardsConfig;
}

/// GENERATOR_X =
/// 5299619240641551281634865583518297030282874472190772894086521144482721001553
pub const GENERATOR_X: CtFq = ct_fq(taceo_ark_babyjubjub::GENERATOR_X.0);

/// GENERATOR_Y =
/// 16950150798460657717958625567821834550301663161624707787222815936182638968203
pub const GENERATOR_Y: CtFq = ct_fq(taceo_ark_babyjubjub::GENERATOR_Y.0);

// Constant-time conditional swap/select helpers for the Montgomery ladder,
// ported from taceo and specialized to the constant-time projective point.

#[inline(always)]
fn conditional_swap(a: &mut CtProjective, b: &mut CtProjective, c: bool) {
    let mask = (c as u64).wrapping_neg(); // all 1s if c, else all 0s
    conditionally_swap_bigint(&mut a.x.0, &mut b.x.0, mask);
    conditionally_swap_bigint(&mut a.y.0, &mut b.y.0, mask);
    conditionally_swap_bigint(&mut a.z.0, &mut b.z.0, mask);
    conditionally_swap_bigint(&mut a.t.0, &mut b.t.0, mask);
}

#[inline(always)]
fn conditional_select(a: &mut CtProjective, b: &CtProjective, c: bool) {
    let mask = (c as u64).wrapping_neg(); // all 1s if c, else all 0s
    conditionally_select_bigint(&mut a.x.0, b.x.0, mask);
    conditionally_select_bigint(&mut a.y.0, b.y.0, mask);
    conditionally_select_bigint(&mut a.z.0, b.z.0, mask);
    conditionally_select_bigint(&mut a.t.0, b.t.0, mask);
}

#[inline(always)]
fn conditionally_select_bigint<const N: usize>(a: &mut BigInt<N>, b: BigInt<N>, mask: u64) {
    for (ai, bi) in a.0.iter_mut().zip(b.0.iter()) {
        *ai ^= mask & (*ai ^ *bi);
    }
}

#[inline(always)]
fn conditionally_swap_bigint<const N: usize>(a: &mut BigInt<N>, b: &mut BigInt<N>, mask: u64) {
    for (ai, bi) in a.0.iter_mut().zip(b.0.iter_mut()) {
        let swap = mask & (*ai ^ *bi);
        *ai ^= swap;
        *bi ^= swap;
    }
}

#[cfg(test)]
mod tests {
    use super::{CtAffine, CtEdwardsConfig, CtProjective};
    use ark_ec::{AdditiveGroup, CurveGroup, twisted_edwards::TECurveConfig};
    use ark_ff::PrimeField;
    use taceo_ark_babyjubjub::{EdwardsConfig as Taceo, EdwardsProjective as TaceoProj};

    fn next_u64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn ct_generator() -> CtProjective {
        CtProjective::from(<CtEdwardsConfig as TECurveConfig>::GENERATOR)
    }

    fn taceo_generator() -> TaceoProj {
        TaceoProj::from(<Taceo as TECurveConfig>::GENERATOR)
    }

    /// Assert two projective points are equal by comparing their canonical
    /// affine coordinates across the two backends.
    fn assert_eq_points(ct: &CtProjective, tc: &TaceoProj) {
        let ca = ct.into_affine();
        let ta = tc.into_affine();
        assert_eq!(ca.x.into_bigint(), ta.x.into_bigint(), "x mismatch");
        assert_eq!(ca.y.into_bigint(), ta.y.into_bigint(), "y mismatch");
    }

    #[test]
    fn generator_matches_taceo() {
        assert_eq_points(&ct_generator(), &taceo_generator());
    }

    #[test]
    fn double_matches_taceo() {
        let mut p_ct = ct_generator();
        let mut p_tc = taceo_generator();
        for _ in 0..8 {
            p_ct.double_in_place();
            p_tc.double_in_place();
            assert_eq_points(&p_ct, &p_tc);
        }
    }

    #[test]
    fn scalar_mul_matches_taceo() {
        // The scalar in `mul_projective` is a raw integer (big-endian bit
        // iteration), so both backends must agree for *any* limbs; no reduction
        // is needed. Cover the edge scalars plus many random ones.
        let g_ct = ct_generator();
        let g_tc = taceo_generator();

        for s in [vec![0u64; 4], {
            let mut v = vec![0u64; 4];
            v[0] = 1;
            v
        }] {
            let pc = <CtEdwardsConfig as TECurveConfig>::mul_projective(&g_ct, &s);
            let pt = <Taceo as TECurveConfig>::mul_projective(&g_tc, &s);
            assert_eq_points(&pc, &pt);
        }

        let mut st = 0x5eed_1234_abcd_0001u64;
        for _ in 0..64 {
            let s = [
                next_u64(&mut st),
                next_u64(&mut st),
                next_u64(&mut st),
                next_u64(&mut st),
            ];
            let pc = <CtEdwardsConfig as TECurveConfig>::mul_projective(&g_ct, &s);
            let pt = <Taceo as TECurveConfig>::mul_projective(&g_tc, &s);
            assert_eq_points(&pc, &pt);
        }
    }

    #[test]
    fn add_matches_taceo() {
        // Build two independent points [a]G and [b]G in each backend, add them,
        // and check the sums agree.
        let g_ct = ct_generator();
        let g_tc = taceo_generator();
        let mut st = 0xfeed_face_0bad_c0deu64;
        for _ in 0..32 {
            let a = [next_u64(&mut st), next_u64(&mut st)];
            let b = [next_u64(&mut st), next_u64(&mut st)];
            let ac = <CtEdwardsConfig as TECurveConfig>::mul_projective(&g_ct, &a);
            let bc = <CtEdwardsConfig as TECurveConfig>::mul_projective(&g_ct, &b);
            let at = <Taceo as TECurveConfig>::mul_projective(&g_tc, &a);
            let bt = <Taceo as TECurveConfig>::mul_projective(&g_tc, &b);
            assert_eq_points(&(ac + bc), &(at + bt));
        }
    }

    #[test]
    fn affine_roundtrip_matches_taceo() {
        // `into_affine` performs a field inversion (constant-time on the CT
        // backend); confirm the affine result matches taceo for many points.
        let g_ct = ct_generator();
        let g_tc = taceo_generator();
        let mut st = 0x0123_4567_89ab_cdefu64;
        for _ in 0..32 {
            let s = [next_u64(&mut st), next_u64(&mut st), next_u64(&mut st)];
            let pc: CtAffine =
                <CtEdwardsConfig as TECurveConfig>::mul_projective(&g_ct, &s).into_affine();
            let pt = <Taceo as TECurveConfig>::mul_projective(&g_tc, &s).into_affine();
            assert_eq!(pc.x.into_bigint(), pt.x.into_bigint());
            assert_eq!(pc.y.into_bigint(), pt.y.into_bigint());
        }
    }
}
