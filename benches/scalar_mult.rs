use babyjubjub_ec::{ProjectivePoint, Scalar};
use criterion::{criterion_group, criterion_main, Criterion};

fn bench_scalar_mult_256bit(c: &mut Criterion) {
    c.bench_function("scalar_mult_256bit", |b| {
        let mut bytes = [0u8; 32];
        bytes.copy_from_slice(&[
            0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB,
            0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67,
            0x89, 0xAB, 0xCD, 0xEF,
        ]);
        let scalar = Scalar::from_bytes_le(&bytes).unwrap();

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR * scalar;
        });
    });
}

fn bench_scalar_mult_large(c: &mut Criterion) {
    c.bench_function("scalar_mult_large", |b| {
        let mut bytes = [0u8; 32];
        // r - 1 (curve order minus 1) - a large scalar
        bytes.copy_from_slice(&[
            0xF0, 0xD3, 0x5C, 0xAE, 0x61, 0x3C, 0xF9, 0x72, 0x9E, 0x67, 0x72, 0x0A, 0xE0, 0x93,
            0x3E, 0x8B, 0x2D, 0xDB, 0x3E, 0xAA, 0x95, 0x16, 0x4C, 0x30, 0xD0, 0xB6, 0x08, 0x70,
            0x3C, 0xE9, 0x65, 0x05,
        ]);
        let scalar = Scalar::from_bytes_le(&bytes).unwrap();

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR * scalar;
        });
    });
}

fn bench_scalar_mult_small(c: &mut Criterion) {
    c.bench_function("scalar_mult_small", |b| {
        let scalar: Scalar = 100u64.into();

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR * scalar;
        });
    });
}

criterion_group!(
    benches,
    bench_scalar_mult_256bit,
    bench_scalar_mult_large,
    bench_scalar_mult_small
);
criterion_main!(benches);
