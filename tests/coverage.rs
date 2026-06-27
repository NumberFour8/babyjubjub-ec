use babyjubjub_ec::{AffinePoint, ProjectivePoint, Scalar};
use group::ff::Field;
use group::Group;
use subtle::ConstantTimeEq;
use rand::SeedableRng;
use rand::rngs::StdRng;

#[test]
fn test_uncovered_randomization() {
    let mut rng = StdRng::seed_from_u64(42);
    // Cover ProjectivePoint::try_random (631-634)
    // Note: ProjectivePoint doesn't implement TryRng-based random in Group trait (it's inherent or other trait)
    // Actually Group trait has `random`, but `try_random` is inherent in src/lib.rs
    let _ = ProjectivePoint::try_random(&mut rng).unwrap();
    
    // Cover Scalar::random (1221-1223, 1500-1505)
    let _ = <Scalar as Field>::random(&mut rng);
    
    // Cover Scalar::try_random (1225-1227, 1507-1511)
    let _ = <Scalar as Field>::try_random(&mut rng).unwrap();
}

#[test]
fn test_uncovered_identity() {
    // Cover ProjectivePoint::identity (640-642)
    let identity = ProjectivePoint::identity();
    assert!(bool::from(identity.is_identity()));
}

#[test]
fn test_uncovered_field_ops() {
    let s = Scalar::from(2u64);
    // Cover Scalar::square (1229-1231)
    let squared = s.square();
    assert_eq!(squared, Scalar::from(4u64));
    
    // Cover Scalar::double (1233-1235)
    let doubled = s.double();
    assert_eq!(doubled, Scalar::from(4u64));
}

#[test]
fn test_uncovered_iterators() {
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
fn test_uncovered_operator_overloads() {
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
fn test_uncovered_ct_eq() {
    let a = AffinePoint::IDENTITY;
    let b = AffinePoint::from(ProjectivePoint::GENERATOR);
    // Cover AffinePoint::ct_eq (1478-1480)
    let _ = a.ct_eq(&b);
}
