#[cfg(feature = "serde")]
mod serde_tests {
    use babyjubjub_ec::{AffinePoint, Scalar};
    use serde_json;

    #[test]
    fn test_scalar_serde_json() {
        let s = Scalar::ONE;
        let serialized = serde_json::to_string(&s).unwrap();
        // Scalar::ONE in big-endian bytes is 00...01
        // Hex representation should be "0000000000000000000000000000000000000000000000000000000000000001"
        assert_eq!(
            serialized,
            "\"0000000000000000000000000000000000000000000000000000000000000001\""
        );

        let deserialized: Scalar = serde_json::from_str(&serialized).unwrap();
        assert_eq!(s, deserialized);
    }

    #[test]
    fn test_affine_point_serde_json() {
        let p = AffinePoint::GENERATOR;
        let serialized = serde_json::to_string(&p).unwrap();

        let deserialized: AffinePoint = serde_json::from_str(&serialized).unwrap();
        assert_eq!(p, deserialized);
        assert!(deserialized.is_on_curve());
        assert!(deserialized.is_in_prime_order_subgroup());
    }

    #[test]
    fn test_scalar_serde_bincode() {
        let s = Scalar::ONE;
        let serialized = bincode::serialize(&s).unwrap();
        assert_eq!(serialized.len(), 32);

        let deserialized: Scalar = bincode::deserialize(&serialized).unwrap();
        assert_eq!(s, deserialized);
    }

    #[test]
    fn test_affine_point_serde_bincode() {
        let p = AffinePoint::GENERATOR;
        let serialized = bincode::serialize(&p).unwrap();
        // 2 coordinates * 32 bytes = 64 bytes
        assert_eq!(serialized.len(), 64);

        let deserialized: AffinePoint = bincode::deserialize(&serialized).unwrap();
        assert_eq!(p, deserialized);
    }

    #[test]
    fn test_scalar_invalid_deser() {
        // Not a hex string
        let res: Result<Scalar, _> = serde_json::from_str("\"not hex\"");
        assert!(res.is_err());

        // Wrong length
        let res: Result<Scalar, _> = serde_json::from_str("\"0102\"");
        assert!(res.is_err());

        // Out of range (modulus + 1)
        // SCALAR_MODULUS = 2736030358979909402780800718157159386076813972158567259200215660948447373041
        // r = 060c89ce5c263405370a08b6d0302b0bab3eedb83920ee0a677297dc392126f1
        let out_of_range = "060c89ce5c263405370a08b6d0302b0bab3eedb83920ee0a677297dc392126f2";
        let res: Result<Scalar, _> = serde_json::from_str(&format!("\"{}\"", out_of_range));
        assert!(res.is_err());
    }
}
