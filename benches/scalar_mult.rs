use babyjubjub_ec::{ProjectivePoint, Scalar};
use criterion::{Criterion, criterion_group, criterion_main};

fn bench_scalar_mult_256bit(c: &mut Criterion) {
    c.bench_function("scalar_mult_256bit", |b| {
        let mut bytes = [0u8; 32];
        bytes.copy_from_slice(&[
            0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB,
            0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67,
            0x89, 0xAB, 0xCD, 0xEF,
        ]);
        // This 256-bit pattern is >= r, so reduce it modulo r rather than decode
        // canonically (`from_bytes_le` would reject it and `unwrap` would panic).
        let scalar = Scalar::reduce_bytes_le(&bytes);

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR * scalar;
        });
    });
}

fn bench_scalar_mult_large(c: &mut Criterion) {
    c.bench_function("scalar_mult_large", |b| {
        let mut bytes = [0u8; 32];
        // An arbitrary large in-range scalar (< r). (Not r-1; just a wide value.)
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

fn bench_mul_fixed_schedule_256bit(c: &mut Criterion) {
    c.bench_function("mul_fixed_schedule_256bit", |b| {
        let mut bytes = [0u8; 32];
        bytes.copy_from_slice(&[
            0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB,
            0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67, 0x89, 0xAB, 0xCD, 0xEF, 0x01, 0x23, 0x45, 0x67,
            0x89, 0xAB, 0xCD, 0xEF,
        ]);
        let scalar = Scalar::reduce_bytes_le(&bytes);

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR.mul_fixed_schedule(&scalar);
        });
    });
}

fn bench_mul_fixed_schedule_large(c: &mut Criterion) {
    c.bench_function("mul_fixed_schedule_large", |b| {
        let mut bytes = [0u8; 32];
        // An arbitrary large in-range scalar (< r). (Not r-1; just a wide value.)
        bytes.copy_from_slice(&[
            0xF0, 0xD3, 0x5C, 0xAE, 0x61, 0x3C, 0xF9, 0x72, 0x9E, 0x67, 0x72, 0x0A, 0xE0, 0x93,
            0x3E, 0x8B, 0x2D, 0xDB, 0x3E, 0xAA, 0x95, 0x16, 0x4C, 0x30, 0xD0, 0xB6, 0x08, 0x70,
            0x3C, 0xE9, 0x65, 0x05,
        ]);
        let scalar = Scalar::from_bytes_le(&bytes).unwrap();

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR.mul_fixed_schedule(&scalar);
        });
    });
}

fn bench_mul_fixed_schedule_small(c: &mut Criterion) {
    c.bench_function("mul_fixed_schedule_small", |b| {
        let scalar: Scalar = 100u64.into();

        b.iter(|| {
            let _ = ProjectivePoint::GENERATOR.mul_fixed_schedule(&scalar);
        });
    });
}

criterion_group!(
    benches,
    bench_scalar_mult_256bit,
    bench_scalar_mult_large,
    bench_scalar_mult_small,
    bench_mul_fixed_schedule_256bit,
    bench_mul_fixed_schedule_large,
    bench_mul_fixed_schedule_small
);
criterion_main!(benches);
