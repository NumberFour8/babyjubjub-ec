//! Public API tests for babyjubjub-ec
//!
//! These tests exercise only the public API and can be run as integration tests.

use babyjubjub_ec::group::ff::{Field, PrimeField};
use babyjubjub_ec::group::{Group, GroupEncoding};
use babyjubjub_ec::subtle::{ConditionallySelectable, ConstantTimeEq};
use babyjubjub_ec::{
    AffinePoint, BabyJubJub, BackendBaseField, GroupRepr, ProjectivePoint, Scalar,
};

// ==================== Identity Tests ====================

#[test]
fn test_affine_point_identity() {
    let identity = AffinePoint::IDENTITY;
    assert!(identity.is_identity());
}

#[test]
fn test_projective_point_identity() {
    let identity = ProjectivePoint::IDENTITY;
    assert!(identity.is_identity());
}

#[test]
fn test_projective_to_affine_identity() {
    let projective = ProjectivePoint::IDENTITY;
    let affine = projective.to_affine();
    assert_eq!(affine, AffinePoint::IDENTITY);
}

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

// ==================== Point Operations Tests ====================

#[test]
fn test_add_identity_projective() {
    let point = ProjectivePoint::GENERATOR;
    let result = point + ProjectivePoint::IDENTITY;
    assert_eq!(result, point);
}

#[test]
fn test_add_refs() {
    let a = ProjectivePoint::GENERATOR;
    let b = ProjectivePoint::IDENTITY;
    let a_ref = &a;
    let b_ref = &b;
    let result = a_ref + b_ref;
    assert_eq!(result, a);
}

#[test]
fn test_point_ops_with_refs() {
    let a = ProjectivePoint::GENERATOR;
    let b = ProjectivePoint::IDENTITY;
    let a_ref = &a;
    let b_ref = &b;
    let result1 = a_ref + b_ref;
    assert_eq!(result1, a);
}

#[test]
fn test_projective_point_neg_ref() {
    let point = ProjectivePoint::GENERATOR;
    let neg_point = -&point;
    let result = point + neg_point;
    assert!(result.is_identity());
}

#[test]
fn test_projective_point_sum() {
    let points = [
        ProjectivePoint::GENERATOR,
        ProjectivePoint::GENERATOR,
        ProjectivePoint::GENERATOR,
    ];
    let sum: ProjectivePoint = points.into_iter().sum();
    let g = ProjectivePoint::GENERATOR;
    let expected = g + g + g;
    assert_eq!(sum, expected);
}

// ==================== Scalar Multiplication Tests ====================

#[test]
fn test_scalar_mult_zero() {
    let scalar_zero: Scalar = 0u64.into();
    let result = ProjectivePoint::GENERATOR * scalar_zero;
    assert_eq!(result, ProjectivePoint::IDENTITY);
}

#[test]
fn test_scalar_mult_one() {
    let scalar_one: Scalar = 1u64.into();
    let result = ProjectivePoint::GENERATOR * scalar_one;
    assert_eq!(result, ProjectivePoint::GENERATOR);
}

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

// ==================== Scalar Field Tests ====================

#[test]
fn test_scalar_add() {
    let a: Scalar = 5u64.into();
    let b: Scalar = 7u64.into();
    let result = a + b;
    let expected: Scalar = 12u64.into();
    assert_eq!(result, expected);
}

#[test]
fn test_scalar_mult() {
    let a: Scalar = 3u64.into();
    let b: Scalar = 4u64.into();
    let result = a * b;
    let expected: Scalar = 12u64.into();
    assert_eq!(result, expected);
}

#[test]
fn test_scalar_neg() {
    let a: Scalar = 5u64.into();
    let neg_a = -a;
    let zero: Scalar = 0u64.into();
    assert_eq!(a + neg_a, zero);
}

#[test]
fn test_scalar_from_u64() {
    let from_u64: Scalar = 1u64.into();
    assert_eq!(from_u64, Scalar::ONE);
}

#[test]
fn test_scalar_from_bytes_be() {
    let mut bytes = [0u8; 32];
    bytes[31] = 42;
    let scalar = Scalar::from_bytes(&bytes).unwrap();
    let expected: Scalar = 42u64.into();
    assert_eq!(scalar, expected);
}

#[test]
fn test_scalar_from_bytes_le_input() {
    let mut bytes = [0u8; 32];
    bytes[0] = 42;
    let scalar = Scalar::from_bytes_le(&bytes).unwrap();
    let expected: Scalar = 42u64.into();
    assert_eq!(scalar, expected);
}

#[test]
fn test_scalar_to_bytes() {
    let scalar: Scalar = 42u64.into();
    let bytes = scalar.to_bytes();
    assert_eq!(bytes[31], 42);
}

#[test]
fn test_scalar_invert() {
    let scalar: Scalar = 5u64.into();
    let inverted = scalar.invert().unwrap();
    let result = scalar * inverted;
    assert_eq!(result, Scalar::ONE);
}

#[test]
fn test_scalar_invert_zero() {
    let scalar = Scalar::ZERO;
    let result = scalar.invert();
    assert!(bool::from(result.is_none()));
}

#[test]
fn test_scalar_square() {
    let scalar: Scalar = 5u64.into();
    let squared = scalar.square();
    let expected: Scalar = 25u64.into();
    assert_eq!(squared, expected);
}

#[test]
fn test_scalar_double() {
    let scalar: Scalar = 5u64.into();
    let doubled = scalar.double();
    let expected: Scalar = 10u64.into();
    assert_eq!(doubled, expected);
}

#[test]
fn test_scalar_is_zero() {
    assert!(Scalar::ZERO.is_zero());
    assert!(!Scalar::ONE.is_zero());
    let non_zero: Scalar = 42u64.into();
    assert!(!non_zero.is_zero());
}

#[test]
fn test_scalar_is_one() {
    assert!(Scalar::ONE.is_one());
    assert!(!Scalar::ZERO.is_one());
    let not_one: Scalar = 42u64.into();
    assert!(!not_one.is_one());
}

#[test]
fn test_scalar_sum() {
    let scalars = [Scalar::from(1u64), Scalar::from(2u64), Scalar::from(3u64)];
    let sum: Scalar = scalars.into_iter().sum();
    let expected = Scalar::from(6u64);
    assert_eq!(sum, expected);
}

#[test]
fn test_scalar_product() {
    let scalars = [Scalar::from(2u64), Scalar::from(3u64), Scalar::from(5u64)];
    let product: Scalar = scalars.into_iter().product();
    let expected = Scalar::from(30u64);
    assert_eq!(product, expected);
}

#[test]
fn test_scalar_sum_empty() {
    let scalars: [Scalar; 0] = [];
    let sum: Scalar = scalars.into_iter().sum();
    assert_eq!(sum, Scalar::ZERO);
}

#[test]
fn test_scalar_sub_neg() {
    let a: Scalar = 42u64.into();
    let zero: Scalar = 0u64.into();
    let result = zero - a;
    assert_eq!(result, -a);
}

#[test]
fn test_scalar_from_bytes_le_known() {
    let mut bytes = [0u8; 32];
    bytes[0] = 0x01;
    let scalar = Scalar::from_bytes_le(&bytes).unwrap();
    assert_eq!(scalar, Scalar::ONE);
}

#[test]
fn test_scalar_is_odd() {
    let even: Scalar = 42u64.into();
    let odd: Scalar = 43u64.into();
    assert!(!bool::from(even.is_odd()));
    assert!(bool::from(odd.is_odd()));
}

#[test]
fn test_scalar_sqrt() {
    let val: Scalar = 16u64.into();
    let sqrt = val.sqrt();
    assert!(bool::from(sqrt.is_some()));
    let result = sqrt.unwrap();
    let squared = result.square();
    assert_eq!(squared, val);
}

#[test]
fn test_scalar_sqrt_two() {
    let val: Scalar = 2u64.into();
    let sqrt = val.sqrt();
    assert!(bool::from(sqrt.is_some()));
    let result = sqrt.unwrap();
    let squared = result.square();
    assert_eq!(squared, val);
}

#[test]
fn test_scalar_sqrt_ratio() {
    // 16/4 = 4 is a square
    let num: Scalar = 16u64.into();
    let den: Scalar = 4u64.into();
    let (is_square, root) = Scalar::sqrt_ratio(&num, &den);
    assert!(bool::from(is_square));
    let ratio = num * den.invert().unwrap();
    assert_eq!(root.square(), ratio);

    // ROOT_OF_UNITY is a non-square
    let (is_square_ns, _) = Scalar::sqrt_ratio(&Scalar::ROOT_OF_UNITY, &Scalar::ONE);
    assert!(!bool::from(is_square_ns));

    // num == 0 => (1, 0)
    let (is_square_n0, root_n0) = Scalar::sqrt_ratio(&Scalar::ZERO, &den);
    assert!(bool::from(is_square_n0));
    assert_eq!(root_n0, Scalar::ZERO);
}

#[test]
fn test_scalar_conditional_select() {
    let a: Scalar = 42u64.into();
    let b: Scalar = 84u64.into();
    let result = Scalar::conditional_select(&a, &b, 0.into());
    assert_eq!(result, a);
    let result = Scalar::conditional_select(&a, &b, 1.into());
    assert_eq!(result, b);
}

#[test]
fn test_ct_eq() {
    let a: Scalar = 42u64.into();
    let b: Scalar = 42u64.into();
    let c: Scalar = 43u64.into();
    assert!(bool::from(a.ct_eq(&b)));
    assert!(!bool::from(a.ct_eq(&c)));
}

#[test]
fn test_scalar_to_bytes_le() {
    let scalar: Scalar = 42u64.into();
    let bytes = scalar.to_bytes_le();
    assert_eq!(bytes[0], 42);
    let reconstructed = Scalar::from_bytes_le(&bytes).unwrap();
    assert_eq!(reconstructed, scalar);
}

#[test]
fn test_scalar_to_bytes_le_zero() {
    let bytes = Scalar::ZERO.to_bytes_le();
    assert_eq!(bytes, [0u8; 32]);
}

#[test]
fn test_scalar_to_bytes_le_one() {
    let bytes = Scalar::ONE.to_bytes_le();
    assert_eq!(bytes[0], 1);
    for b in bytes.iter().skip(1) {
        assert_eq!(*b, 0);
    }
}

#[test]
fn test_two_inv() {
    let result = Scalar::TWO_INV * Scalar::from(2u64);
    assert_eq!(result, Scalar::ONE);
}

// ==================== Affine/Projective Conversion Tests ====================

#[test]
fn test_affine_to_projective() {
    let affine = AffinePoint::GENERATOR;
    let projective = ProjectivePoint::from(affine);
    assert_eq!(projective.to_affine(), affine);
}

#[test]
fn test_affine_projective_round_trip() {
    let affine = AffinePoint::GENERATOR;
    let projective: ProjectivePoint = affine.into();
    let affine_back: AffinePoint = projective.into();
    assert_eq!(affine, affine_back);
}

#[test]
fn test_projective_to_affine_generator() {
    let projective = ProjectivePoint::GENERATOR;
    let affine = projective.to_affine();
    assert_eq!(affine, AffinePoint::GENERATOR);
}

// ==================== AffinePoint Tests ====================

#[test]
fn test_affine_point_x() {
    let affine = AffinePoint::GENERATOR;
    assert_eq!(affine.x(), affine.x);
}

#[test]
fn test_affine_point_y() {
    let affine = AffinePoint::GENERATOR;
    assert_eq!(affine.y(), affine.y);
}

#[test]
fn test_affine_point_is_identity() {
    assert!(AffinePoint::IDENTITY.is_identity());
    assert!(!AffinePoint::GENERATOR.is_identity());
}

#[test]
fn test_affine_point_x_is_odd() {
    let affine = AffinePoint::GENERATOR;
    let _ = affine.x_is_odd();
}

#[test]
fn test_affine_point_conditional_select() {
    let a = AffinePoint::IDENTITY;
    let b = AffinePoint::GENERATOR;
    let result = AffinePoint::conditional_select(&a, &b, 0.into());
    assert_eq!(result, a);
    let result = AffinePoint::conditional_select(&a, &b, 1.into());
    assert_eq!(result, b);
}

// ==================== GroupEncoding Tests ====================

#[test]
fn test_group_encoding_to_bytes() {
    let point = ProjectivePoint::GENERATOR;
    let repr = point.to_bytes();
    assert_eq!(repr.as_ref().len(), 32);
}

#[test]
fn test_group_encoding_from_bytes() {
    // All zeros is invalid (not in prime-order subgroup)
    let bytes = GroupRepr([0u8; 32]);
    let result = ProjectivePoint::from_bytes(&bytes);
    assert!(bool::from(result.is_none()));
}

#[test]
fn test_group_encoding_round_trip_identity() {
    let identity = ProjectivePoint::IDENTITY;
    let bytes = identity.to_bytes();
    assert_eq!(bytes.as_ref().len(), 32);

    let decoded = ProjectivePoint::from_bytes(&bytes);
    assert!(bool::from(decoded.is_some()));
    let decoded_point = decoded.unwrap();
    assert!(decoded_point.is_identity());
}

#[test]
fn test_group_encoding_round_trip_generator() {
    let generator = ProjectivePoint::GENERATOR;
    let bytes = generator.to_bytes();

    let decoded = ProjectivePoint::from_bytes(&bytes);
    assert!(bool::from(decoded.is_some()));
    let decoded_point = decoded.unwrap();

    assert_eq!(decoded_point.to_affine(), generator.to_affine());
}

#[test]
fn test_group_encoding_from_bytes_unchecked() {
    // The all-zero encoding has y = 0, which decodes to a cofactor point of
    // order 4: on the curve but NOT in the prime-order subgroup. This is the
    // case that actually distinguishes the two decoders -- `from_bytes` runs the
    // subgroup check and must reject it, while `from_bytes_unchecked` skips that
    // check and must accept it.
    let bytes = GroupRepr([0u8; 32]);

    // Checked decoding rejects the small-subgroup point.
    assert!(bool::from(ProjectivePoint::from_bytes(&bytes).is_none()));

    // Unchecked decoding accepts it, and the recovered point is indeed outside
    // the prime-order subgroup -- so the two paths really do behave differently.
    let decoded = ProjectivePoint::from_bytes_unchecked(&bytes);
    assert!(bool::from(decoded.is_some()));
    assert!(!decoded.unwrap().is_in_prime_order_subgroup());
}

#[test]
fn test_group_encoding_from_bytes_valid_point() {
    let scalar: Scalar = 42u64.into();
    let point = ProjectivePoint::GENERATOR * scalar;
    let bytes = point.to_bytes();

    let decoded = ProjectivePoint::from_bytes(&bytes);
    assert!(bool::from(decoded.is_some()));

    let decoded_affine = decoded.unwrap().to_affine();
    let point_affine = point.to_affine();
    assert_eq!(decoded_affine, point_affine);
}

#[test]
fn test_group_encoding_from_bytes_sign_bit() {
    // A point and its negation share the same y-coordinate and differ only in
    // the sign of x, which the compressed encoding packs into the top bit of
    // the final byte. Decoding must honor that bit and recover the correct
    // point rather than its negation.
    let g = ProjectivePoint::GENERATOR;
    let neg_g = -g;

    let g_bytes = g.to_bytes();
    let neg_bytes = neg_g.to_bytes();

    // Identical y => only the packed x-sign bit (bit 7 of the last byte) differs.
    assert_eq!(g_bytes.as_ref()[..31], neg_bytes.as_ref()[..31]);
    assert_eq!(g_bytes.as_ref()[31] ^ neg_bytes.as_ref()[31], 0x80);

    // Each encoding decodes back to its own point, with the sign preserved.
    let g_decoded = ProjectivePoint::from_bytes(&g_bytes).unwrap();
    let neg_decoded = ProjectivePoint::from_bytes(&neg_bytes).unwrap();
    assert_eq!(g_decoded.to_affine(), g.to_affine());
    assert_eq!(neg_decoded.to_affine(), neg_g.to_affine());
    assert_ne!(g_decoded.to_affine(), neg_decoded.to_affine());
}

#[test]
fn test_group_encoding_from_bytes_invalid_point() {
    // Should not panic
    let mut bytes = [0u8; 32];
    bytes[0] = 2;
    let _ = ProjectivePoint::from_bytes(&GroupRepr(bytes));
}

#[test]
fn test_point_encoding_is_canonical() {
    let g = ProjectivePoint::GENERATOR;
    let bytes = g.to_bytes();
    assert_eq!(bytes.as_ref().len(), 32);
    assert!(bool::from(ProjectivePoint::from_bytes(&bytes).is_some()));

    // Non-canonical spare bit must be rejected
    let mut mutated = bytes;
    mutated.0[31] |= 0x40;
    assert!(bool::from(ProjectivePoint::from_bytes(&mutated).is_none()));
}

// ==================== GroupRepr Tests ====================

#[test]
fn test_group_repr_default() {
    let repr = GroupRepr::default();
    assert_eq!(repr.as_ref(), &[0u8; 32]);
}

#[test]
fn test_group_repr_as_ref() {
    let repr = GroupRepr([42u8; 32]);
    let bytes: &[u8] = repr.as_ref();
    assert_eq!(bytes.len(), 32);
    assert_eq!(bytes[0], 42);
}

#[test]
fn test_group_repr_as_mut() {
    let mut repr = GroupRepr([0u8; 32]);
    let bytes: &mut [u8] = repr.as_mut();
    bytes[0] = 42;
    assert_eq!(repr.as_ref()[0], 42);
}

// ==================== Validation Tests ====================

#[test]
fn test_on_curve_and_subgroup_helpers() {
    // Generator is on-curve and in the prime-order subgroup
    assert!(AffinePoint::GENERATOR.is_on_curve());
    assert!(AffinePoint::GENERATOR.is_in_prime_order_subgroup());
    assert!(ProjectivePoint::GENERATOR.is_in_prime_order_subgroup());
    assert!(ProjectivePoint::GENERATOR.is_on_curve());

    // Valid point via new (checked constructor)
    let g = AffinePoint::GENERATOR;
    assert!(AffinePoint::new(g.x, g.y).is_some());
}

// ==================== Conversion Tests ====================

#[test]
fn test_from_affine_ref_to_projective() {
    let affine = AffinePoint::GENERATOR;
    let projective: ProjectivePoint = (&affine).into();
    assert_eq!(projective.to_affine(), affine);
}

#[test]
fn test_from_projective_ref_to_affine() {
    let projective = ProjectivePoint::GENERATOR;
    let affine: AffinePoint = (&projective).into();
    assert_eq!(affine, projective.to_affine());
}

// ==================== Scalar Constants Tests ====================

#[test]
fn test_scalar_two_inv() {
    let result = Scalar::TWO_INV * Scalar::from(2u64);
    assert_eq!(result, Scalar::ONE);
}

#[test]
fn test_scalar_num_bits() {
    assert_eq!(Scalar::NUM_BITS, 251);
}

#[test]
fn test_scalar_capacity() {
    assert_eq!(Scalar::CAPACITY, 250);
}

#[test]
fn test_scalar_s() {
    assert_eq!(Scalar::S, 4);
}

// ==================== Cofactor Clearing Tests ====================

#[test]
fn test_mul_with_cofactor_clear_matches_eight_s_on_subgroup_point() {
    // On a prime-order point, mul_with_cofactor_clear(s) == [8s]P for any scalar,
    // including large ones where the integer 8*s wraps past r. Test both a z==1
    // base (GENERATOR) and a z!=1 base (G+G+G) so the to_backend_unchecked
    // extended-coordinate path is exercised for non-unit z.
    let g = ProjectivePoint::GENERATOR;
    let g3 = g + g + g; // in the subgroup, but z != 1
    let eight = Scalar::from(8u64);
    for base in [g, g3] {
        for s in [
            Scalar::from(1u64),
            Scalar::from(7u64),
            Scalar::from(123456u64),
            Scalar::from(u64::MAX) * Scalar::from(0x9e37_79b9_7f4a_7c15u64),
            Scalar::ZERO - Scalar::ONE, // r - 1: 8*s wraps
        ] {
            let cleared = base.mul_with_cofactor_clear(&s);
            assert_eq!(
                cleared.to_affine(),
                (base * (s * eight)).to_affine(),
                "mul_with_cofactor_clear(s) must equal [8s]P"
            );
            assert!(
                cleared.is_in_prime_order_subgroup(),
                "result must be in the prime-order subgroup"
            );
        }
    }
}

// ==================== Constant-Time Equality Tests ====================

#[test]
fn test_ct_eq_projective_scaled_and_identity() {
    let g = ProjectivePoint::GENERATOR;
    // A scaled (non-unit-z) representation of the same point must compare equal.
    let p = g + g; // z != 1
    let scaled = (p + g) - g; // same affine point as p, un-normalized
    assert!(bool::from(p.ct_eq(&scaled)));
    assert!(!bool::from(p.ct_eq(&g)));
    // Identity in two representations.
    assert!(bool::from(ProjectivePoint::IDENTITY.ct_eq(&(g + (-g)))));
}

#[test]
fn test_ct_eq_projective_does_not_panic_on_zero_z() {
    // A z==0 point is reachable via struct instantiation / Default. ct_eq must not panic
    // (the cross-multiplication implementation needs no inversion).
    let zero = BackendBaseField::from(0u64);
    let one = BackendBaseField::from(1u64);
    let invalid = ProjectivePoint {
        x: zero,
        y: one,
        z: zero,
    };
    // Comparing equal invalid points returns true; comparing to the generator false.
    assert!(bool::from(invalid.ct_eq(&invalid)));
    assert!(!bool::from(invalid.ct_eq(&ProjectivePoint::GENERATOR)));
    // Default is also z==0 and must not panic.
    let _ = ProjectivePoint::default().ct_eq(&ProjectivePoint::default());
}

// ==================== Encoding Canonicality / Malleability Tests ====================

#[test]
fn test_from_bytes_rejects_noncanonical_identity_sign_bit() {
    // The identity (x==0) has its x-sign flag in bit 7 of the last byte. Since
    // x == -x == 0, an encoding with that bit set is a second, non-canonical
    // encoding of the identity and must be rejected (non-malleability).
    let canonical = ProjectivePoint::IDENTITY.to_bytes();
    assert_eq!(
        canonical.as_ref()[31] & 0x80,
        0,
        "canonical identity has sign bit 0"
    );
    // The canonical identity still decodes.
    assert!(bool::from(
        ProjectivePoint::from_bytes(&canonical).is_some()
    ));

    let mut malleated = canonical;
    malleated.0[31] |= 0x80;
    assert!(
        bool::from(ProjectivePoint::from_bytes(&malleated).is_none()),
        "identity with x-sign bit set must be rejected by from_bytes"
    );
    assert!(
        bool::from(ProjectivePoint::from_bytes_unchecked(&malleated).is_none()),
        "identity with x-sign bit set must be rejected by from_bytes_unchecked too"
    );
}

// ==================== is_on_curve (projective, z != 1) Tests ====================

#[test]
fn test_is_on_curve_scaled_projective() {
    // A genuine scaled projective point (z != 1) on the curve returns true.
    let p = ProjectivePoint::GENERATOR + ProjectivePoint::GENERATOR;
    assert!(p.is_on_curve());
    // Perturbing a coordinate takes it off the curve (exercises the false branch
    // of the projective curve equation with z != 1).
    let bad = ProjectivePoint {
        x: p.x + BackendBaseField::from(1u64),
        y: p.y,
        z: p.z,
    };
    assert!(!bad.is_on_curve());
}

// ==================== Random Point Tests ====================

#[test]
fn test_group_random_is_in_prime_order_subgroup() {
    use rand::{SeedableRng, rngs::StdRng};
    let mut rng = StdRng::seed_from_u64(7);
    let mut seen_distinct = false;
    let first = <ProjectivePoint as Group>::random(&mut rng);
    for _ in 0..16 {
        let p = <ProjectivePoint as Group>::random(&mut rng);
        assert!(
            p.is_in_prime_order_subgroup(),
            "Group::random must yield a prime-order-subgroup point"
        );
        if !bool::from(p.ct_eq(&first)) {
            seen_distinct = true;
        }
    }
    assert!(seen_distinct, "random must produce distinct points");
}

// ==================== PrimeField Repr Round-Trip Tests ====================

#[test]
fn test_prime_field_repr_round_trip() {
    for s in [
        Scalar::ZERO,
        Scalar::ONE,
        Scalar::from(42u64),
        Scalar::ZERO - Scalar::ONE,
    ] {
        let repr = s.to_repr();
        let back = Scalar::from_repr(repr);
        assert!(bool::from(back.is_some()));
        assert_eq!(back.unwrap(), s);
    }
    // from_repr rejects a non-canonical (>= r) encoding.
    assert!(bool::from(Scalar::from_repr([0xFFu8; 32].into()).is_none()));
}

// ==================== Assign-Operator Tests ====================

#[test]
fn test_scalar_assign_operators() {
    let a = Scalar::from(20u64);
    let b = Scalar::from(7u64);

    let mut x = a;
    x += b;
    assert_eq!(x, a + b);
    let mut x = a;
    x += &b;
    assert_eq!(x, a + b);

    let mut x = a;
    x -= b;
    assert_eq!(x, a - b);
    let mut x = a;
    x -= &b;
    assert_eq!(x, a - b);

    let mut x = a;
    x *= b;
    assert_eq!(x, a * b);
    let mut x = a;
    x *= &b;
    assert_eq!(x, a * b);
}

#[test]
fn test_projective_assign_operators() {
    let g = ProjectivePoint::GENERATOR;
    let s = Scalar::from(5u64);

    let mut p = g;
    p += g;
    assert_eq!(p, g + g);
    let mut p = g;
    p += &g;
    assert_eq!(p, g + g);

    let mut p = g + g;
    p -= g;
    assert_eq!(p, g);
    let mut p = g + g;
    p -= &g;
    assert_eq!(p, g);

    let mut p = g;
    p *= s;
    assert_eq!(p, g * s);
    let mut p = g;
    p *= &s;
    assert_eq!(p, g * s);
}

// ==================== PrimeCurve Marker Test ====================

#[test]
fn test_babyjubjub_implements_prime_curve() {
    // Compile-time check that BabyJubJub satisfies the PrimeCurve bound, plus a
    // sanity check that its ORDER constant is exposed via Curve.
    fn assert_prime_curve<C: elliptic_curve::PrimeCurve>() {}
    assert_prime_curve::<BabyJubJub>();
    let _ = <BabyJubJub as elliptic_curve::Curve>::ORDER;
}

// ==================== CofactorGroup Tests ====================

#[test]
fn test_projective_point_implements_cofactor_group() {
    use babyjubjub_ec::group::cofactor::CofactorGroup;
    use babyjubjub_ec::group::prime::PrimeGroup;

    // Compile-time check that ProjectivePoint satisfies the CofactorGroup and
    // PrimeGroup bounds, with the subgroup type being ProjectivePoint itself.
    fn assert_cofactor_group<G>()
    where
        G: CofactorGroup<Subgroup = ProjectivePoint> + PrimeGroup,
    {
    }
    assert_cofactor_group::<ProjectivePoint>();

    // Functional sanity check through the public trait surface: the generator is
    // torsion free, `into_subgroup` accepts it, and clearing the cofactor keeps
    // it inside the prime-order subgroup.
    let g = ProjectivePoint::GENERATOR;
    assert!(bool::from(g.is_torsion_free()));
    assert!(bool::from(g.into_subgroup().is_some()));
    assert!(g.clear_cofactor().is_in_prime_order_subgroup());
}

// ==================== CurveArithmetic Tests ====================

#[test]
fn test_babyjubjub_implements_curve_arithmetic() {
    // Compile-time check that BabyJubJub satisfies the full CurveArithmetic bound.
    fn assert_curve_arithmetic<C: elliptic_curve::CurveArithmetic>() {}
    assert_curve_arithmetic::<BabyJubJub>();

    // The associated types are exactly the crate's wrapper types (this only
    // type-checks if they are identical).
    let _scalar: <BabyJubJub as elliptic_curve::CurveArithmetic>::Scalar = Scalar::ONE;
    let _affine: <BabyJubJub as elliptic_curve::CurveArithmetic>::AffinePoint =
        AffinePoint::GENERATOR;
    let _projective: <BabyJubJub as elliptic_curve::CurveArithmetic>::ProjectivePoint =
        ProjectivePoint::GENERATOR;
}

#[test]
fn test_scalar_u256_round_trip() {
    use elliptic_curve::bigint::U256;
    use elliptic_curve::bigint::modular::Retrieve;
    use elliptic_curve::scalar::FromUintUnchecked;

    for s in [
        Scalar::ZERO,
        Scalar::ONE,
        Scalar::from(42u64),
        Scalar::ZERO - Scalar::ONE, // r - 1, the largest canonical scalar
    ] {
        let u: U256 = s.into();
        // `retrieve` agrees with the `Into<U256>` conversion.
        assert_eq!(s.retrieve(), u);
        // `from_uint_unchecked` round-trips canonical (< r) integers.
        assert_eq!(Scalar::from_uint_unchecked(u), s);
    }
}

#[test]
fn test_scalar_scalarprimitive_round_trip() {
    use elliptic_curve::ScalarPrimitive;

    for s in [Scalar::ZERO, Scalar::ONE, Scalar::from(123u64)] {
        let sv: ScalarPrimitive<BabyJubJub> = s.into();
        let back: Scalar = sv.into();
        assert_eq!(back, s);
    }
}

#[test]
fn test_scalar_fieldbytes_round_trip() {
    use elliptic_curve::FieldBytes;

    for s in [
        Scalar::ZERO,
        Scalar::ONE,
        Scalar::from(99u64),
        Scalar::ZERO - Scalar::ONE,
    ] {
        let fb: FieldBytes<BabyJubJub> = s.into();
        // `Into<FieldBytes>` matches `PrimeField::to_repr`.
        assert_eq!(fb, s.to_repr());
        assert_eq!(Scalar::from_repr(fb).unwrap(), s);
    }
}

#[test]
fn test_nonzero_scalar_conversions() {
    use elliptic_curve::NonZeroScalar;

    // Zero is rejected, both via `new` and `TryFrom` (with `elliptic_curve::Error`).
    assert!(bool::from(
        NonZeroScalar::<BabyJubJub>::new(Scalar::ZERO).is_none()
    ));
    assert!(NonZeroScalar::<BabyJubJub>::try_from(Scalar::ZERO).is_err());

    // A non-zero scalar is accepted and round-trips losslessly.
    let s = Scalar::from(7u64);
    let nz = NonZeroScalar::<BabyJubJub>::try_from(s).expect("non-zero scalar");
    let back: Scalar = nz.into();
    assert_eq!(back, s);
}

#[test]
fn test_scalar_reduce_u256_and_fieldbytes() {
    use elliptic_curve::FieldBytes;
    use elliptic_curve::bigint::U256;
    use elliptic_curve::ops::Reduce;

    // Values < r are returned unchanged.
    assert_eq!(
        <Scalar as Reduce<U256>>::reduce(U256::from(5u64)),
        Scalar::from(5u64)
    );

    // A value >= r is reduced; the result is canonical, hence idempotent.
    let reduced = <Scalar as Reduce<U256>>::reduce(U256::MAX);
    let reduced_u: U256 = reduced.into();
    assert_eq!(<Scalar as Reduce<U256>>::reduce(reduced_u), reduced);

    // `reduce_bytes` reduces a non-canonical little-endian encoding too.
    let fb: FieldBytes<BabyJubJub> = [0xFFu8; 32].into();
    let r1 = <Scalar as Reduce<U256>>::reduce_bytes(&fb);
    let r1_bytes: FieldBytes<BabyJubJub> = r1.into();
    assert_eq!(<Scalar as Reduce<U256>>::reduce_bytes(&r1_bytes), r1);
}

#[test]
fn test_scalar_is_high() {
    use elliptic_curve::scalar::IsHigh;

    assert!(!bool::from(Scalar::ZERO.is_high()));
    assert!(!bool::from(Scalar::ONE.is_high()));
    // r - 1 is the largest scalar and is unambiguously in the upper half.
    assert!(bool::from((Scalar::ZERO - Scalar::ONE).is_high()));
}

#[test]
fn test_scalar_shr_assign() {
    // `ShrAssign<usize>` is required by the 0.13.8 `CurveArithmetic::Scalar`
    // bound: a right shift is floor-division of the canonical integer value by
    // `2^n`.
    let mut s = Scalar::from(180u64);
    s >>= 2usize;
    assert_eq!(s, Scalar::from(45u64)); // 180 / 4

    // Shifting by 0 bits is a no-op.
    let mut t = Scalar::from(12_345u64);
    t >>= 0usize;
    assert_eq!(t, Scalar::from(12_345u64));

    // Property: for an even value 2k, (2k) >> 1 == k.
    let k = Scalar::from(7u64);
    let mut two_k = k + k;
    two_k >>= 1usize;
    assert_eq!(two_k, k);
}

#[test]
#[allow(clippy::op_ref)] // the point of this test is to exercise the reference operator impls
fn test_scalar_times_point_equals_point_times_scalar() {
    let s = Scalar::from(123_456_789u64);
    let g = ProjectivePoint::GENERATOR;
    let expected = g * s;

    // Scalar * ProjectivePoint, in every operand/reference combination.
    assert_eq!(s * g, expected);
    assert_eq!(&s * g, expected);
    assert_eq!(s * &g, expected);
    assert_eq!(&s * &g, expected);

    // Scalar * AffinePoint matches, in both operand orders.
    let ga = AffinePoint::GENERATOR;
    assert_eq!(s * ga, expected);
    assert_eq!(ga * s, expected);
}

#[test]
fn test_linear_combination_two_terms() {
    use elliptic_curve::ops::LinearCombination;

    let g = ProjectivePoint::GENERATOR;
    let q = g + g; // 2G
    let a = Scalar::from(3u64);
    let b = Scalar::from(5u64);

    // 0.13.8's `LinearCombination::lincomb` takes the two point/scalar pairs as
    // separate arguments and computes `x * k + y * l`.
    let combined = ProjectivePoint::lincomb(&g, &a, &q, &b);
    assert_eq!(combined, g * a + q * b);
}

#[test]
fn test_curve_group_to_affine_and_batch_normalize() {
    use elliptic_curve::BatchNormalize;

    let g = ProjectivePoint::GENERATOR;
    let p = g + g; // z != 1

    // `group::Curve::to_affine` matches the inherent to_affine.
    assert_eq!(
        <ProjectivePoint as elliptic_curve::group::Curve>::to_affine(&p),
        p.to_affine()
    );

    // BatchNormalize of an array equals per-element normalization.
    let arr = [g, p, g + p];
    let normalized =
        <ProjectivePoint as BatchNormalize<[ProjectivePoint; 3]>>::batch_normalize(&arr);
    for (point, affine) in arr.iter().zip(normalized.iter()) {
        assert_eq!(*affine, point.to_affine());
    }

    // 0.13.8 has no `CurveAffine`; identity / generator / to-curve are exposed
    // via inherent methods, the `GENERATOR`/`IDENTITY` constants, and `From`.
    let affine = p.to_affine();
    assert_eq!(ProjectivePoint::from(affine), p);
    assert!(AffinePoint::IDENTITY.is_identity());
    assert!(!AffinePoint::GENERATOR.is_identity());
}

#[test]
fn test_affine_coordinates() {
    use elliptic_curve::FieldBytes;
    use elliptic_curve::point::AffineCoordinates;

    let g = AffinePoint::GENERATOR;
    // `AffineCoordinates::x` returns the x-coordinate as `FieldBytes`; 0.13.8's
    // trait has no `y`, so the inherent `coordinates` getter supplies both
    // coordinates as `FieldBytes` (the inherent `x`/`y` return field elements).
    let gx = AffineCoordinates::x(&g);
    let (gx2, gy) = g.coordinates();
    assert_eq!(gx, gx2);
    let rebuilt = AffinePoint::from_coordinates(&gx, &gy);
    assert!(bool::from(rebuilt.is_some()));
    assert_eq!(rebuilt.unwrap(), g);

    // (0, 0) is not on the curve and must be rejected.
    let zero = FieldBytes::<BabyJubJub>::default();
    assert!(bool::from(
        AffinePoint::from_coordinates(&zero, &zero).is_none()
    ));
}

#[test]
fn test_nonidentity_conversions() {
    use elliptic_curve::point::NonIdentity;

    // The identity is rejected for both point types.
    assert!(NonIdentity::<ProjectivePoint>::try_from(ProjectivePoint::IDENTITY).is_err());
    assert!(NonIdentity::<AffinePoint>::try_from(AffinePoint::IDENTITY).is_err());

    // A non-identity projective point is accepted and round-trips.
    let g = ProjectivePoint::GENERATOR;
    let nz = NonIdentity::<ProjectivePoint>::try_from(g).expect("generator is non-identity");
    let back: ProjectivePoint = nz.into();
    assert_eq!(back, g);

    // Same for affine.
    let ga = AffinePoint::GENERATOR;
    let nza = NonIdentity::<AffinePoint>::try_from(ga).expect("generator is non-identity");
    let backa: AffinePoint = nza.into();
    assert_eq!(backa, ga);
}

#[test]
fn test_affine_point_default_is_identity() {
    // `Default` for both point types is now the group identity (required so the
    // generic `NonIdentity` machinery treats it as the identity sentinel).
    assert_eq!(AffinePoint::default(), AffinePoint::IDENTITY);
    assert_eq!(ProjectivePoint::default(), ProjectivePoint::IDENTITY);
}
