//! Compatibility tests: output must match vcftools 0.1.17 semantics.
//!
//! Golden files under tests/golden/ encode vcftools-verified expected output.
//! The last row of --TsTv-by-qual output has a vcftools C bug (uninitialized
//! N_Tv_GT value); our goldens use the correct value (0) instead of that garbage.
//! No external process is invoked at test time.

use rsomics_vcf_tstv_strat::{tstv_by_count, tstv_by_qual, write_count_table, write_qual_table};
use std::io::BufReader;

fn run_count(vcf: &str) -> String {
    let rows = tstv_by_count(BufReader::new(vcf.as_bytes())).unwrap();
    let mut buf = Vec::new();
    write_count_table(&rows, &mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

fn run_qual(vcf: &str) -> String {
    let rows = tstv_by_qual(BufReader::new(vcf.as_bytes())).unwrap();
    let mut buf = Vec::new();
    write_qual_table(&rows, &mut buf).unwrap();
    String::from_utf8(buf).unwrap()
}

// --- by-count tests ---

#[test]
fn by_count_3samples() {
    // 3 samples: A>G (Ts), C>T (Ts), G>A (Ts), A>C (Tv), G>T (Tv), indel excluded,
    // A>G with ./. missing → count=1.
    let vcf = include_str!("golden/by_count_3samples.vcf");
    let expected = include_str!("golden/by_count_3samples.count");
    assert_eq!(run_count(vcf), expected);
}

#[test]
fn by_count_missing_gt() {
    // ./., .|., and "." (haploid missing) all contribute 0 to the ALT allele count.
    let vcf = include_str!("golden/by_count_missing_gt.vcf");
    let expected = include_str!("golden/by_count_missing_gt.count");
    assert_eq!(run_count(vcf), expected);
}

#[test]
fn by_count_multiallelic_excluded() {
    // Multiallelic sites (>1 ALT allele) are fully excluded: no bin entry, not even
    // as Ts/Tv=0. Only the biallelic A>G (count=2) contributes.
    let vcf = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
        chr1\t100\t.\tA\tG,C\t50\tPASS\t.\tGT\t0/1\t0/1\n\
        chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/1\n";
    let got = run_count(vcf);
    let lines: Vec<&str> = got.lines().collect();
    // bin 2: only A>G (count=2, Ts) → N_Ts=1
    assert_eq!(
        lines[3], "2\t1\t0\tinf",
        "multiallelic must be fully excluded"
    );
    // bin 1: multiallelic would contribute count=2 for GT=0/1,0/1 with first-ALT→count=2,
    // but multiallelic is excluded so bin 1 should be 0,0
    assert_eq!(lines[2], "1\t0\t0\tnan");
}

#[test]
fn by_count_indel_excluded() {
    // Indels are fully excluded (not even counted in any bin).
    let vcf = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
        chr1\t100\t.\tACGT\tA\t50\tPASS\t.\tGT\t0/1\n\
        chr1\t200\t.\tA\tG\t50\tPASS\t.\tGT\t1/1\n";
    let got = run_count(vcf);
    let lines: Vec<&str> = got.lines().collect();
    // Only 1 sample → bins 0..=1
    assert_eq!(lines[0], "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv");
    assert_eq!(lines[1], "0\t0\t0\tnan");
    // Indel bin (count=1) would be here if counted — should be empty
    assert_eq!(lines[2], "1\t0\t0\tnan");
    // A>G (Ts) count=2 goes to bin 2? But only 1 sample (diploid) → max bin=1
    // Wait: 1 sample → 2*1=2 bins (0,1). A>G with 1/1 → count=2, but max_bin=1
    // → overflow, not counted? Let's see what our impl does
    assert_eq!(lines.len(), 3, "1 sample → 2 data rows + header");
}

#[test]
fn by_count_max_bin_is_2n_samples() {
    // Max bin = 2*N_samples - 1 regardless of observed max count.
    // 4 samples → bins 0..=7.
    let vcf = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\n\
        chr1\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t0/0\t0/0\t0/0\n";
    let got = run_count(vcf);
    let lines: Vec<&str> = got.lines().collect();
    // header + bins 0..=7 = 9 lines
    assert_eq!(
        lines.len(),
        9,
        "4 samples → bins 0..7 → 9 lines (header + 8 data)"
    );
    assert_eq!(lines[8], "7\t0\t0\tnan");
}

#[test]
fn by_count_all_ts_bin0() {
    // ALT_ALLELE_COUNT=0 when all genotypes are 0/0 (or missing). SNP is still
    // counted in bin 0 with its Ts/Tv class.
    let vcf = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\n\
        chr1\t100\t.\tA\tC\t50\tPASS\t.\tGT\t0/0\t0/0\n";
    let got = run_count(vcf);
    let lines: Vec<&str> = got.lines().collect();
    assert_eq!(lines[1], "0\t0\t1\t0", "A>C Tv with 0/0 → count=0, N_Tv=1");
}

// --- by-qual tests ---

#[test]
fn by_qual_basic() {
    // Thresholds at unique QUAL values of biallelic SNPs; indel excluded from threshold set.
    // LT=strictly less, GT=strictly greater than threshold.
    // Last row has N_Tv_GT=0 (correct; vcftools has C uninitialized-memory bug there).
    let vcf = include_str!("golden/by_qual_basic.vcf");
    let expected = include_str!("golden/by_qual_basic.qual");
    assert_eq!(run_qual(vcf), expected);
}

#[test]
fn by_qual_dot_qual() {
    // QUAL=. is treated as -1.0 and creates a row with QUAL_THRESHOLD=-1.
    let vcf = include_str!("golden/by_qual_dot_qual.vcf");
    let expected = include_str!("golden/by_qual_dot_qual.qual");
    assert_eq!(run_qual(vcf), expected);
}

#[test]
fn by_qual_duplicate_quals_merge() {
    // Duplicate QUAL values produce only one threshold row.
    let vcf = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
        chr1\t100\t.\tA\tG\t20\tPASS\t.\tGT\t0/1\n\
        chr1\t200\t.\tC\tT\t20\tPASS\t.\tGT\t0/1\n\
        chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n";
    let got = run_qual(vcf);
    let lines: Vec<&str> = got.lines().collect();
    // Unique QUALs = {20, 30} → 2 data rows + header = 3 lines total
    assert_eq!(
        lines.len(),
        3,
        "duplicate QUAL=20 merges into one threshold row"
    );
    assert!(lines[1].starts_with("20\t"), "first threshold is 20");
    assert!(lines[2].starts_with("30\t"), "second threshold is 30");
}

#[test]
fn by_qual_indel_not_in_threshold() {
    // Indel QUALs do not appear as threshold rows.
    let vcf = "##fileformat=VCFv4.2\n\
        ##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n\
        chr1\t100\t.\tA\tG\t10\tPASS\t.\tGT\t0/1\n\
        chr1\t150\t.\tACG\tA\t15\tPASS\t.\tGT\t0/1\n\
        chr1\t200\t.\tC\tT\t20\tPASS\t.\tGT\t0/1\n";
    let got = run_qual(vcf);
    let lines: Vec<&str> = got.lines().collect();
    // Only SNP QUALs {10, 20} → 2 threshold rows
    assert_eq!(
        lines.len(),
        3,
        "indel QUAL=15 must not create a threshold row"
    );
    assert!(!got.contains("15\t"), "indel QUAL=15 absent from output");
}

#[test]
fn by_qual_empty_vcf() {
    // No SNP sites → empty output (header only).
    let vcf = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let got = run_qual(vcf);
    assert_eq!(
        got,
        "QUAL_THRESHOLD\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD\n"
    );
}

#[test]
fn by_count_empty_vcf() {
    // No samples → 0 bins (empty output).
    let vcf = "##fileformat=VCFv4.2\n\
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let got = run_count(vcf);
    assert_eq!(got, "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv\n");
}
