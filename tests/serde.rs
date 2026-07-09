#[cfg(feature = "serde")]
mod tests {
    use babyjubjub_ec::{Scalar, AffinePoint, ProjectivePoint, GroupRepr};
    use group::GroupEncoding;

    #[test]
    fn test_scalar_serde() {
        let scalar = Scalar::ONE;
        let serialized = serde_json::to_string(&scalar).unwrap();
        let deserialized: Scalar = serde_json::from_str(&serialized).unwrap();
        assert_eq!(scalar, deserialized);
    }

    #[test]
    fn test_affine_point_serde() {
        let point = AffinePoint::GENERATOR;
        let serialized = serde_json::to_string(&point).unwrap();
        let deserialized: AffinePoint = serde_json::from_str(&serialized).unwrap();
        assert_eq!(point, deserialized);
    }

    #[test]
    fn test_group_repr_serde() {
        let point = ProjectivePoint::GENERATOR;
        let repr = point.to_bytes();
        let serialized = serde_json::to_string(&repr).unwrap();
        let deserialized: GroupRepr = serde_json::from_str(&serialized).unwrap();
        assert_eq!(repr, deserialized);
    }
}
