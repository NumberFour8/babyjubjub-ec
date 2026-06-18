//! Public API tests for babyjubjub-ec
//!
//! These tests exercise only the public API and can be run as integration tests.

use babyjubjub_ec::{
    AffinePoint, GroupRepr, ProjectivePoint, Scalar,
};
use group::{Group, GroupEncoding};
use group::ff::{Field, PrimeField};
use subtle::{ConditionallySelectable, ConstantTimeEq};

// ==================== Identity Tests ====================

#[test]
fn test_affine_point_identity() {
    let identity = AffinePoint::IDENTITY;
    assert!(bool::from(identity.is_identity()));
}

#[test]
fn test_projective_point_identity() {
    let identity = ProjectivePoint::IDENTITY;
    assert!(bool::from(identity.is_identity()));
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
    let result = &a + &b;
    assert_eq!(result, a);
}

#[test]
fn test_point_ops_with_refs() {
    let a = ProjectivePoint::GENERATOR;
    let b = ProjectivePoint::IDENTITY;
    let result1 = &a + &b;
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
    let points = vec![
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
    assert_eq!(g.mul_fixed_schedule(&Scalar::ONE).to_affine(), g.to_affine());
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
    let scalars = vec![Scalar::from(1u64), Scalar::from(2u64), Scalar::from(3u64)];
    let sum: Scalar = scalars.into_iter().sum();
    let expected = Scalar::from(6u64);
    assert_eq!(sum, expected);
}

#[test]
fn test_scalar_product() {
    let scalars = vec![Scalar::from(2u64), Scalar::from(3u64), Scalar::from(5u64)];
    let product: Scalar = scalars.into_iter().product();
    let expected = Scalar::from(30u64);
    assert_eq!(product, expected);
}

#[test]
fn test_scalar_sum_empty() {
    let scalars: Vec<Scalar> = vec![];
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
    let gen = ProjectivePoint::GENERATOR;
    let bytes = gen.to_bytes();

    let decoded = ProjectivePoint::from_bytes(&bytes);
    assert!(bool::from(decoded.is_some()));
    let decoded_point = decoded.unwrap();

    assert_eq!(decoded_point.to_affine(), gen.to_affine());
}

#[test]
fn test_group_encoding_from_bytes_unchecked() {
    let gen = ProjectivePoint::GENERATOR;
    let bytes = gen.to_bytes();
    let result = ProjectivePoint::from_bytes_unchecked(&bytes);
    assert!(bool::from(result.is_some()));
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
    let gen = ProjectivePoint::GENERATOR;
    let bytes = gen.to_bytes();

    let decoded = ProjectivePoint::from_bytes(&bytes);
    assert!(bool::from(decoded.is_some()));

    let decoded_affine = decoded.unwrap().to_affine();
    let gen_affine = gen.to_affine();
    assert_eq!(decoded_affine, gen_affine);
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
    assert!(bool::from(AffinePoint::GENERATOR.is_on_curve()));
    assert!(bool::from(AffinePoint::GENERATOR.is_in_prime_order_subgroup()));
    assert!(bool::from(ProjectivePoint::GENERATOR.is_in_prime_order_subgroup()));
    assert!(bool::from(ProjectivePoint::GENERATOR.is_on_curve()));

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