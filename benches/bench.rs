use criterion::{Criterion, criterion_group, criterion_main};
use rsomics_vcf_tstv_strat::{tstv_by_count, tstv_by_qual};
use std::io::BufReader;

fn make_large_vcf(n_variants: usize, n_samples: usize) -> String {
    let snp_pairs = [
        ("A", "G"),
        ("C", "T"),
        ("G", "A"),
        ("A", "C"),
        ("G", "T"),
        ("A", "T"),
    ];
    let mut buf = String::from("##fileformat=VCFv4.2\n");
    buf.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    buf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for s in 0..n_samples {
        buf.push_str(&format!("\ts{s}"));
    }
    buf.push('\n');

    for i in 0..n_variants {
        let chrom = format!("chr{}", (i % 22) + 1);
        let pos = i * 100 + 1;
        let qual = 10 + (i % 90);
        let (r, a) = snp_pairs[i % snp_pairs.len()];
        buf.push_str(&format!("{chrom}\t{pos}\t.\t{r}\t{a}\t{qual}\tPASS\t.\tGT"));
        for s in 0..n_samples {
            let gt = match (i + s) % 3 {
                0 => "0/0",
                1 => "0/1",
                _ => "1/1",
            };
            buf.push('\t');
            buf.push_str(gt);
        }
        buf.push('\n');
    }
    buf
}

fn bench_by_count(c: &mut Criterion) {
    let vcf = make_large_vcf(100_000, 10);
    c.bench_function("tstv_by_count_100k_10samples", |b| {
        b.iter(|| {
            let rows = tstv_by_count(BufReader::new(vcf.as_bytes())).unwrap();
            std::hint::black_box(rows);
        });
    });
}

fn bench_by_qual(c: &mut Criterion) {
    let vcf = make_large_vcf(100_000, 10);
    c.bench_function("tstv_by_qual_100k_10samples", |b| {
        b.iter(|| {
            let rows = tstv_by_qual(BufReader::new(vcf.as_bytes())).unwrap();
            std::hint::black_box(rows);
        });
    });
}

criterion_group!(benches, bench_by_count, bench_by_qual);
criterion_main!(benches);
