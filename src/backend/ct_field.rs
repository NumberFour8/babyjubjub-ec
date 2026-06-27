//! Constant-time `ark_ff::FpConfig<4>` field backends bridged to the vendored,
//! formally-verified fiat-crypto arithmetic. Compiled only with `--features fiat`.
//!
//! Both backends expose the field as a standard `ark_ff::Fp<_, 4>` whose every
//! arithmetic operation delegates to branch-free fiat-crypto code with masked
//! (data-independent) reduction, and whose inversion is the constant-time
//! Bernstein-Yang `divstep` chain. Because `ark_ff::Fp` stores its limbs in
//! Montgomery form with `R = 2^256` — exactly fiat-crypto's `word-by-word-
//! montgomery` domain — the bridge is a direct limb reinterpretation with no
//! conversion. All field constants are sourced (by copying the identical
//! Montgomery limbs) from the existing arkworks/taceo backend, so they are
//! correct by construction and additionally checked for equality in tests.
//!
//! - [`CtFq`]: BabyJubJub base field (the BN254 scalar prime); fixes the M1
//!   field-level timing leak in curve add/double/mul.
//! - [`CtFr`]: BabyJubJub scalar field; provides constant-time inversion (M2).

#![allow(clippy::too_many_arguments)]

use ark_ff::{BigInt, Fp, FpConfig, MontBackend, SqrtPrecomputation};
use core::marker::PhantomData;

use super::fiat::{bjjfr, bn254fq};

/// taceo's (variable-time) configs for the same two primes. Used *only* as a
/// compile-time source of field constants: identical prime and identical
/// Montgomery radix (`R = 2^256`) mean the limb representation is bit-identical,
/// so reinterpreting their constant limbs yields the correct constants for the
/// constant-time configs.
type TaceoFqCfg = MontBackend<taceo_ark_babyjubjub::FqConfig, 4>;
type TaceoFrCfg = MontBackend<taceo_ark_babyjubjub::FrConfig, 4>;

/// Raw-limb bridge over one field's vendored fiat-crypto functions.
///
/// Every method works on plain `[u64; 4]`/`[u64; 5]` limb arrays (Montgomery
/// form for field elements), hiding fiat's per-field domain newtypes. This lets
/// the arithmetic and the Bernstein-Yang inversion driver be written once,
/// generically over the field.
trait FiatBackend {
    /// Bernstein-Yang `divstep` iteration count for this prime:
    /// `(49 * bits + 57) / 17` with `bits = floor(log2(m)) + 1`.
    const ITERATIONS: usize;

    fn mul(out: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]);
    fn add(out: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]);
    fn sub(out: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]);
    fn square(out: &mut [u64; 4], a: &[u64; 4]);
    fn opp(out: &mut [u64; 4], a: &[u64; 4]);
    /// Montgomery domain -> standard (integer) domain.
    fn from_montgomery(out: &mut [u64; 4], a: &[u64; 4]);
    /// Standard (integer) domain -> Montgomery domain.
    fn to_montgomery(out: &mut [u64; 4], a: &[u64; 4]);
    /// Montgomery representation of `1`.
    fn set_one(out: &mut [u64; 4]);
    /// The modulus in saturated (5-limb, two's-complement) form.
    fn msat(out: &mut [u64; 5]);
    /// Precomputed `R^2 / 2^(ITERATIONS)` correction factor (Montgomery form).
    fn divstep_precomp(out: &mut [u64; 4]);
    /// Branch-free select: `out = if c == 0 { a } else { b }`.
    fn selectznz(out: &mut [u64; 4], c: u8, a: &[u64; 4], b: &[u64; 4]);
    /// Non-zero indicator: returns `0` iff `a` is the zero element.
    fn nonzero(a: &[u64; 4]) -> u64;
    /// One Bernstein-Yang division step (see fiat-crypto's `divstep`).
    fn divstep(
        out1: &mut u64,
        out2: &mut [u64; 5],
        out3: &mut [u64; 5],
        out4: &mut [u64; 4],
        out5: &mut [u64; 4],
        d: u64,
        f: &[u64; 5],
        g: &[u64; 5],
        v: &[u64; 4],
        r: &[u64; 4],
    );
}

/// Generates a [`FiatBackend`] impl for a marker type by wrapping the field's
/// vendored fiat-crypto functions (handling fiat's domain newtypes).
macro_rules! fiat_backend {
    (
        $backend:ident, $module:ident, $iterations:expr,
        $mont:ident, $nonmont:ident,
        $mul:ident, $add:ident, $sub:ident, $square:ident, $opp:ident,
        $from_mont:ident, $to_mont:ident, $set_one:ident, $msat:ident,
        $divstep:ident, $divstep_precomp:ident, $selectznz:ident, $nonzero:ident
    ) => {
        struct $backend;

        impl FiatBackend for $backend {
            const ITERATIONS: usize = $iterations;

            #[inline]
            fn mul(out: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$mul(&mut o, &$module::$mont(*a), &$module::$mont(*b));
                *out = o.0;
            }
            #[inline]
            fn add(out: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$add(&mut o, &$module::$mont(*a), &$module::$mont(*b));
                *out = o.0;
            }
            #[inline]
            fn sub(out: &mut [u64; 4], a: &[u64; 4], b: &[u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$sub(&mut o, &$module::$mont(*a), &$module::$mont(*b));
                *out = o.0;
            }
            #[inline]
            fn square(out: &mut [u64; 4], a: &[u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$square(&mut o, &$module::$mont(*a));
                *out = o.0;
            }
            #[inline]
            fn opp(out: &mut [u64; 4], a: &[u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$opp(&mut o, &$module::$mont(*a));
                *out = o.0;
            }
            #[inline]
            fn from_montgomery(out: &mut [u64; 4], a: &[u64; 4]) {
                let mut o = $module::$nonmont([0u64; 4]);
                $module::$from_mont(&mut o, &$module::$mont(*a));
                *out = o.0;
            }
            #[inline]
            fn to_montgomery(out: &mut [u64; 4], a: &[u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$to_mont(&mut o, &$module::$nonmont(*a));
                *out = o.0;
            }
            #[inline]
            fn set_one(out: &mut [u64; 4]) {
                let mut o = $module::$mont([0u64; 4]);
                $module::$set_one(&mut o);
                *out = o.0;
            }
            #[inline]
            fn msat(out: &mut [u64; 5]) {
                $module::$msat(out);
            }
            #[inline]
            fn divstep_precomp(out: &mut [u64; 4]) {
                $module::$divstep_precomp(out);
            }
            #[inline]
            fn selectznz(out: &mut [u64; 4], c: u8, a: &[u64; 4], b: &[u64; 4]) {
                $module::$selectznz(out, c, a, b);
            }
            #[inline]
            fn nonzero(a: &[u64; 4]) -> u64 {
                let mut o = 0u64;
                $module::$nonzero(&mut o, a);
                o
            }
            #[inline]
            fn divstep(
                out1: &mut u64,
                out2: &mut [u64; 5],
                out3: &mut [u64; 5],
                out4: &mut [u64; 4],
                out5: &mut [u64; 4],
                d: u64,
                f: &[u64; 5],
                g: &[u64; 5],
                v: &[u64; 4],
                r: &[u64; 4],
            ) {
                $module::$divstep(out1, out2, out3, out4, out5, d, f, g, v, r);
            }
        }
    };
}

fiat_backend!(
    FqBackend,
    bn254fq,
    735,
    fiat_bn254fq_montgomery_domain_field_element,
    fiat_bn254fq_non_montgomery_domain_field_element,
    fiat_bn254fq_mul,
    fiat_bn254fq_add,
    fiat_bn254fq_sub,
    fiat_bn254fq_square,
    fiat_bn254fq_opp,
    fiat_bn254fq_from_montgomery,
    fiat_bn254fq_to_montgomery,
    fiat_bn254fq_set_one,
    fiat_bn254fq_msat,
    fiat_bn254fq_divstep,
    fiat_bn254fq_divstep_precomp,
    fiat_bn254fq_selectznz,
    fiat_bn254fq_nonzero
);

fiat_backend!(
    FrBackend,
    bjjfr,
    726,
    fiat_bjjfr_montgomery_domain_field_element,
    fiat_bjjfr_non_montgomery_domain_field_element,
    fiat_bjjfr_mul,
    fiat_bjjfr_add,
    fiat_bjjfr_sub,
    fiat_bjjfr_square,
    fiat_bjjfr_opp,
    fiat_bjjfr_from_montgomery,
    fiat_bjjfr_to_montgomery,
    fiat_bjjfr_set_one,
    fiat_bjjfr_msat,
    fiat_bjjfr_divstep,
    fiat_bjjfr_divstep_precomp,
    fiat_bjjfr_selectznz,
    fiat_bjjfr_nonzero
);

/// Constant-time modular inversion via fiat-crypto's Bernstein-Yang `divstep`
/// chain. `a` is a Montgomery-form element; returns `a^-1` (Montgomery form),
/// or `None` iff `a == 0`. This is a faithful port of fiat-crypto's reference
/// `inversion/c/inversion_template.c` driver, which runs a fixed number of
/// `divstep`s and is therefore branch-free in the value of `a`.
fn ct_inverse<B: FiatBackend>(a_mont: &[u64; 4]) -> Option<[u64; 4]> {
    if B::nonzero(a_mont) == 0 {
        return None;
    }

    let mut precomp = [0u64; 4];
    B::divstep_precomp(&mut precomp);

    let mut d: u64 = 1;
    let mut f = [0u64; 5];
    B::msat(&mut f);

    // g = a in the standard (non-Montgomery) domain, zero-extended to 5 limbs.
    let mut g = [0u64; 5];
    {
        let mut ns = [0u64; 4];
        B::from_montgomery(&mut ns, a_mont);
        g[..4].copy_from_slice(&ns);
    }

    let mut v = [0u64; 4];
    let mut r = [0u64; 4];
    B::set_one(&mut r);

    // Scratch outputs for the alternating divstep buffers.
    let mut o1 = 0u64;
    let mut o2 = [0u64; 5];
    let mut o3 = [0u64; 5];
    let mut o4 = [0u64; 4];
    let mut o5 = [0u64; 4];

    let iters = B::ITERATIONS;
    let mut i = 0;
    while i < iters - (iters % 2) {
        B::divstep(
            &mut o1, &mut o2, &mut o3, &mut o4, &mut o5, d, &f, &g, &v, &r,
        );
        B::divstep(
            &mut d, &mut f, &mut g, &mut v, &mut r, o1, &o2, &o3, &o4, &o5,
        );
        i += 2;
    }
    if iters % 2 == 1 {
        B::divstep(
            &mut o1, &mut o2, &mut o3, &mut o4, &mut o5, d, &f, &g, &v, &r,
        );
        v = o4;
        f = o2;
    }

    // Sign correction: if the running `f` ended negative (top bit of its top
    // saturated limb), negate `v`; then multiply by the precomputed factor.
    let mut neg = [0u64; 4];
    B::opp(&mut neg, &v);
    let sign = (f[4] >> 63) as u8;
    let mut v_signed = [0u64; 4];
    B::selectznz(&mut v_signed, sign, &v, &neg);

    let mut out = [0u64; 4];
    B::mul(&mut out, &v_signed, &precomp);
    Some(out)
}

/// Generates the constant-time `FpConfig<4>` impl plus the `Fp` type alias for
/// one field, delegating arithmetic to `$backend` and sourcing all field
/// constants from `$taceo` (identical Montgomery limbs).
macro_rules! impl_ct_fp_config {
    ($cfg:ident, $fp:ident, $backend:ident, $taceo:ty) => {
        /// Constant-time field configuration (see module docs).
        pub struct $cfg;

        /// Constant-time prime field element backed by fiat-crypto arithmetic.
        pub type $fp = Fp<$cfg, 4>;

        impl FpConfig<4> for $cfg {
            const MODULUS: BigInt<4> = <$taceo as FpConfig<4>>::MODULUS;
            const GENERATOR: $fp = Fp(<$taceo as FpConfig<4>>::GENERATOR.0, PhantomData);
            const ZERO: $fp = Fp(<$taceo as FpConfig<4>>::ZERO.0, PhantomData);
            const ONE: $fp = Fp(<$taceo as FpConfig<4>>::ONE.0, PhantomData);
            const TWO_ADICITY: u32 = <$taceo as FpConfig<4>>::TWO_ADICITY;
            const TWO_ADIC_ROOT_OF_UNITY: $fp = Fp(
                <$taceo as FpConfig<4>>::TWO_ADIC_ROOT_OF_UNITY.0,
                PhantomData,
            );
            const SMALL_SUBGROUP_BASE: Option<u32> = <$taceo as FpConfig<4>>::SMALL_SUBGROUP_BASE;
            const SMALL_SUBGROUP_BASE_ADICITY: Option<u32> =
                <$taceo as FpConfig<4>>::SMALL_SUBGROUP_BASE_ADICITY;
            const LARGE_SUBGROUP_ROOT_OF_UNITY: Option<$fp> =
                match <$taceo as FpConfig<4>>::LARGE_SUBGROUP_ROOT_OF_UNITY {
                    Some(v) => Some(Fp(v.0, PhantomData)),
                    None => None,
                };
            const SQRT_PRECOMP: Option<SqrtPrecomputation<$fp>> =
                match <$taceo as FpConfig<4>>::SQRT_PRECOMP {
                    Some(SqrtPrecomputation::TonelliShanks {
                        two_adicity,
                        quadratic_nonresidue_to_trace,
                        trace_of_modulus_minus_one_div_two,
                    }) => Some(SqrtPrecomputation::TonelliShanks {
                        two_adicity,
                        quadratic_nonresidue_to_trace: Fp(
                            quadratic_nonresidue_to_trace.0,
                            PhantomData,
                        ),
                        trace_of_modulus_minus_one_div_two,
                    }),
                    Some(SqrtPrecomputation::Case3Mod4 {
                        modulus_plus_one_div_four,
                    }) => Some(SqrtPrecomputation::Case3Mod4 {
                        modulus_plus_one_div_four,
                    }),
                    _ => None,
                };

            #[inline]
            fn add_assign(a: &mut $fp, b: &$fp) {
                let mut o = [0u64; 4];
                <$backend>::add(&mut o, &a.0.0, &b.0.0);
                a.0.0 = o;
            }
            #[inline]
            fn sub_assign(a: &mut $fp, b: &$fp) {
                let mut o = [0u64; 4];
                <$backend>::sub(&mut o, &a.0.0, &b.0.0);
                a.0.0 = o;
            }
            #[inline]
            fn double_in_place(a: &mut $fp) {
                let mut o = [0u64; 4];
                <$backend>::add(&mut o, &a.0.0, &a.0.0);
                a.0.0 = o;
            }
            #[inline]
            fn neg_in_place(a: &mut $fp) {
                let mut o = [0u64; 4];
                <$backend>::opp(&mut o, &a.0.0);
                a.0.0 = o;
            }
            #[inline]
            fn mul_assign(a: &mut $fp, b: &$fp) {
                let mut o = [0u64; 4];
                <$backend>::mul(&mut o, &a.0.0, &b.0.0);
                a.0.0 = o;
            }
            #[inline]
            fn square_in_place(a: &mut $fp) {
                let mut o = [0u64; 4];
                <$backend>::square(&mut o, &a.0.0);
                a.0.0 = o;
            }

            fn sum_of_products<const T: usize>(a: &[$fp; T], b: &[$fp; T]) -> $fp {
                let mut acc = [0u64; 4];
                let mut i = 0;
                while i < T {
                    let mut prod = [0u64; 4];
                    <$backend>::mul(&mut prod, &a[i].0.0, &b[i].0.0);
                    let mut s = [0u64; 4];
                    <$backend>::add(&mut s, &acc, &prod);
                    acc = s;
                    i += 1;
                }
                Fp(BigInt(acc), PhantomData)
            }

            #[inline]
            fn inverse(a: &$fp) -> Option<$fp> {
                ct_inverse::<$backend>(&a.0.0).map(|limbs| Fp(BigInt(limbs), PhantomData))
            }

            fn from_bigint(b: BigInt<4>) -> Option<$fp> {
                if b >= Self::MODULUS {
                    return None;
                }
                let mut o = [0u64; 4];
                <$backend>::to_montgomery(&mut o, &b.0);
                Some(Fp(BigInt(o), PhantomData))
            }

            #[inline]
            fn into_bigint(a: $fp) -> BigInt<4> {
                let mut o = [0u64; 4];
                <$backend>::from_montgomery(&mut o, &a.0.0);
                BigInt(o)
            }
        }
    };
}

impl_ct_fp_config!(CtFqConfig, CtFq, FqBackend, TaceoFqCfg);
impl_ct_fp_config!(CtFrConfig, CtFr, FrBackend, TaceoFrCfg);

/// Reinterprets raw Montgomery limbs (e.g. copied from the taceo backend's
/// curve constants, which share the identical `R = 2^256` representation) as a
/// constant-time base-field element. `const` so it can build curve constants.
pub(crate) const fn ct_fq(limbs: BigInt<4>) -> CtFq {
    Fp(limbs, PhantomData)
}

/// Scalar-field analogue of [`ct_fq`].
pub(crate) const fn ct_fr(limbs: BigInt<4>) -> CtFr {
    Fp(limbs, PhantomData)
}

#[cfg(test)]
mod tests {
    /// Deterministic SplitMix64 PRNG so the equivalence tests need no rng dep
    /// and are fully reproducible.
    fn next_u64(state: &mut u64) -> u64 {
        *state = state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = *state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn rand_bytes(state: &mut u64) -> [u8; 32] {
        let mut b = [0u8; 32];
        let mut k = 0;
        while k < 4 {
            b[k * 8..k * 8 + 8].copy_from_slice(&next_u64(state).to_le_bytes());
            k += 1;
        }
        b
    }

    const N_ITERS: usize = 256;

    /// Generates the cross-backend equivalence tests for one field: every
    /// arithmetic op, inversion, sqrt and the field constants must match the
    /// existing (variable-time) taceo backend exactly.
    macro_rules! equiv_tests {
        ($modname:ident, $ct:ty, $ctcfg:ty, $tc:ty, $tccfg:ty) => {
            mod $modname {
                use super::{N_ITERS, rand_bytes};
                use ark_ff::{AdditiveGroup, Field, FpConfig, PrimeField, Zero};

                type Ct = $ct;
                type Tc = $tc;

                fn ct_from(t: Tc) -> Ct {
                    Ct::from_bigint(t.into_bigint()).unwrap()
                }

                #[test]
                fn arithmetic_matches() {
                    let mut st = 0x1234_5678_9abc_def0u64;
                    for _ in 0..N_ITERS {
                        let ba = rand_bytes(&mut st);
                        let bb = rand_bytes(&mut st);
                        let at = Tc::from_le_bytes_mod_order(&ba);
                        let bt = Tc::from_le_bytes_mod_order(&bb);
                        let ac = ct_from(at);
                        let bc = ct_from(bt);

                        assert_eq!((at + bt).into_bigint(), (ac + bc).into_bigint());
                        assert_eq!((at - bt).into_bigint(), (ac - bc).into_bigint());
                        assert_eq!((bt - at).into_bigint(), (bc - ac).into_bigint());
                        assert_eq!((at * bt).into_bigint(), (ac * bc).into_bigint());
                        assert_eq!(at.square().into_bigint(), ac.square().into_bigint());
                        assert_eq!((-at).into_bigint(), (-ac).into_bigint());
                        assert_eq!(at.double().into_bigint(), ac.double().into_bigint());
                        assert_eq!(
                            at.inverse().map(|x| x.into_bigint()),
                            ac.inverse().map(|x| x.into_bigint())
                        );

                        // `from_le_bytes_mod_order` and bigint round-trips match.
                        assert_eq!(
                            Ct::from_le_bytes_mod_order(&ba).into_bigint(),
                            at.into_bigint()
                        );
                        assert_eq!(
                            Ct::from_bigint(ac.into_bigint()).unwrap().into_bigint(),
                            ac.into_bigint()
                        );
                    }
                }

                #[test]
                fn inverse_edge_cases() {
                    assert!(Ct::ZERO.inverse().is_none());
                    assert_eq!(
                        Ct::ONE.inverse().unwrap().into_bigint(),
                        Ct::ONE.into_bigint()
                    );

                    let mut st = 0x0bad_f00d_dead_beefu64;
                    for _ in 0..N_ITERS {
                        let t = Tc::from_le_bytes_mod_order(&rand_bytes(&mut st));
                        if t.is_zero() {
                            continue;
                        }
                        let c = ct_from(t);
                        assert_eq!(
                            (c.inverse().unwrap() * c).into_bigint(),
                            Ct::ONE.into_bigint()
                        );
                    }
                }

                #[test]
                fn sum_of_products_matches() {
                    let mut st = 0xfeed_face_cafe_d00du64;
                    for _ in 0..N_ITERS {
                        let mut at = [Tc::ZERO; 3];
                        let mut bt = [Tc::ZERO; 3];
                        let mut ac = [Ct::ZERO; 3];
                        let mut bc = [Ct::ZERO; 3];
                        for k in 0..3 {
                            at[k] = Tc::from_le_bytes_mod_order(&rand_bytes(&mut st));
                            bt[k] = Tc::from_le_bytes_mod_order(&rand_bytes(&mut st));
                            ac[k] = ct_from(at[k]);
                            bc[k] = ct_from(bt[k]);
                        }
                        let expected =
                            <$tccfg as FpConfig<4>>::sum_of_products::<3>(&at, &bt).into_bigint();
                        let got =
                            <$ctcfg as FpConfig<4>>::sum_of_products::<3>(&ac, &bc).into_bigint();
                        assert_eq!(expected, got);
                    }
                }

                #[test]
                fn sqrt_matches() {
                    let mut st = 0xabcd_1234_5678_9abcu64;
                    for _ in 0..N_ITERS {
                        let t = Tc::from_le_bytes_mod_order(&rand_bytes(&mut st));
                        let sq_t = t.square();
                        let sq_c = ct_from(t).square();
                        let root_t = sq_t.sqrt().unwrap();
                        let root_c = sq_c.sqrt().unwrap();
                        assert_eq!(root_t.square().into_bigint(), sq_t.into_bigint());
                        assert_eq!(root_c.square().into_bigint(), sq_c.into_bigint());
                    }
                }

                #[test]
                fn constants_match() {
                    assert_eq!(
                        <$ctcfg as FpConfig<4>>::MODULUS,
                        <$tccfg as FpConfig<4>>::MODULUS
                    );
                    assert_eq!(
                        <$ctcfg as FpConfig<4>>::ZERO.into_bigint(),
                        <$tccfg as FpConfig<4>>::ZERO.into_bigint()
                    );
                    assert_eq!(
                        <$ctcfg as FpConfig<4>>::ONE.into_bigint(),
                        <$tccfg as FpConfig<4>>::ONE.into_bigint()
                    );
                    assert_eq!(
                        <$ctcfg as FpConfig<4>>::GENERATOR.into_bigint(),
                        <$tccfg as FpConfig<4>>::GENERATOR.into_bigint()
                    );
                    assert_eq!(
                        <$ctcfg as FpConfig<4>>::TWO_ADICITY,
                        <$tccfg as FpConfig<4>>::TWO_ADICITY
                    );
                    assert_eq!(
                        <$ctcfg as FpConfig<4>>::TWO_ADIC_ROOT_OF_UNITY.into_bigint(),
                        <$tccfg as FpConfig<4>>::TWO_ADIC_ROOT_OF_UNITY.into_bigint()
                    );
                }
            }
        };
    }

    equiv_tests!(
        fq,
        crate::backend::ct_field::CtFq,
        crate::backend::ct_field::CtFqConfig,
        taceo_ark_babyjubjub::Fq,
        crate::backend::ct_field::TaceoFqCfg
    );
    equiv_tests!(
        fr,
        crate::backend::ct_field::CtFr,
        crate::backend::ct_field::CtFrConfig,
        taceo_ark_babyjubjub::Fr,
        crate::backend::ct_field::TaceoFrCfg
    );
}
