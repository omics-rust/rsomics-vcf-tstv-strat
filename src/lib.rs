#![allow(clippy::cast_precision_loss)]

use std::io::BufRead;

use rsomics_common::{Result, RsomicsError};
use serde::Serialize;

/// One row of `--TsTv-by-count` output.
#[derive(Debug, Clone, Serialize)]
pub struct CountRow {
    pub alt_allele_count: u32,
    pub n_ts: u64,
    pub n_tv: u64,
    pub ts_tv: String,
}

/// One row of `--TsTv-by-qual` output.
#[derive(Debug, Clone, Serialize)]
pub struct QualRow {
    pub qual_threshold: String,
    pub n_ts_lt: u64,
    pub n_tv_lt: u64,
    pub ts_tv_lt: String,
    pub n_ts_gt: u64,
    pub n_tv_gt: u64,
    pub ts_tv_gt: String,
}

/// Format Ts/Tv ratio as vcftools 0.1.17 does (C `%g`, precision 6).
///
/// nan when both zero; inf when Tv=0 and Ts>0.
pub fn fmt_tstv(n_ts: u64, n_tv: u64) -> String {
    if n_tv == 0 && n_ts == 0 {
        return "nan".to_owned();
    }
    if n_tv == 0 {
        return "inf".to_owned();
    }
    format_g6(n_ts as f64 / n_tv as f64)
}

/// Reproduce C `printf("%g", x)` with precision 6: fixed notation when the
/// post-rounding decimal exponent is in `-4..6`, scientific otherwise; trailing
/// zeros and a bare `.` are stripped; lowercase nan/inf; 0 prints as "0".
fn format_g6(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_owned();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_owned();
    }
    if x == 0.0 {
        return "0".to_owned();
    }
    const PRECISION: i32 = 6;
    let sci = format!("{:.*e}", (PRECISION - 1) as usize, x);
    let (mantissa, exp_str) = sci.split_once('e').unwrap();
    let exp: i32 = exp_str.parse().unwrap();
    if (-4..PRECISION).contains(&exp) {
        let decimals = (PRECISION - 1 - exp).max(0) as usize;
        strip_g(format!("{x:.decimals$}"))
    } else {
        let sign = if exp < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", strip_g(mantissa.to_owned()), sign, exp.abs())
    }
}

fn strip_g(mut s: String) -> String {
    if s.contains('.') {
        while s.ends_with('0') {
            s.pop();
        }
        if s.ends_with('.') {
            s.pop();
        }
    }
    s
}

/// Classify a biallelic SNP as transition (true) or transversion (false).
///
/// Returns None if REF or ALT is not a single DNA base, or if the pair is not a valid SNP.
fn classify_snp(r#ref: &[u8], alt: &[u8]) -> Option<bool> {
    if r#ref.len() != 1 || alt.len() != 1 {
        return None;
    }
    let r = r#ref[0].to_ascii_uppercase();
    let a = alt[0].to_ascii_uppercase();
    match (r, a) {
        (b'A', b'G') | (b'G', b'A') | (b'C', b'T') | (b'T', b'C') => Some(true),
        (b'A', b'C')
        | (b'C', b'A')
        | (b'A', b'T')
        | (b'T', b'A')
        | (b'G', b'C')
        | (b'C', b'G')
        | (b'G', b'T')
        | (b'T', b'G') => Some(false),
        _ => None,
    }
}

/// Sum ALT allele calls from a GT field string (e.g. "0/1", "1|1", "./.", ".").
///
/// Counts any allele index ≥ 1 as one ALT call. Missing alleles (`.`) contribute 0.
fn count_alt_alleles(gt: &str) -> u32 {
    gt.split(['/', '|'])
        .filter(|a| *a != "." && a.chars().all(|c| c.is_ascii_digit()))
        .filter(|a| a.parse::<u32>().unwrap_or(0) >= 1)
        .count() as u32
}

/// Parse a VCF record line into (ref_bytes, first_alt_bytes, format_idx, sample_fields_start).
///
/// Returns None for header lines. Fails loud on malformed data lines.
fn parse_line(line: &str) -> Result<Option<VcfRecord<'_>>> {
    if line.starts_with('#') {
        return Ok(None);
    }
    let mut fields = line.splitn(10, '\t');
    let _chrom = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing CHROM".to_owned()))?;
    let _pos = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing POS".to_owned()))?;
    let _id = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing ID".to_owned()))?;
    let r#ref = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing REF".to_owned()))?;
    let alt_field = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing ALT".to_owned()))?;
    let qual_str = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing QUAL".to_owned()))?;
    let _filter = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing FILTER".to_owned()))?;
    let _info = fields
        .next()
        .ok_or_else(|| RsomicsError::InvalidInput("missing INFO".to_owned()))?;

    // FORMAT field (may be absent for sites-only VCF)
    let format = fields.next();
    // Remaining: sample columns joined by remaining tabs (we re-split later)
    let samples_tail = fields.next().unwrap_or("");

    // Only biallelic sites: single ALT allele, no comma
    let is_biallelic = !alt_field.contains(',');
    let first_alt = alt_field.split(',').next().unwrap_or(alt_field);

    let qual: Option<f64> = if qual_str == "." {
        Some(-1.0)
    } else {
        qual_str.parse::<f64>().ok()
    };

    Ok(Some(VcfRecord {
        r#ref,
        first_alt,
        is_biallelic,
        qual,
        format,
        samples_tail,
    }))
}

struct VcfRecord<'a> {
    r#ref: &'a str,
    first_alt: &'a str,
    is_biallelic: bool,
    qual: Option<f64>,
    format: Option<&'a str>,
    samples_tail: &'a str,
}

impl<'a> VcfRecord<'a> {
    /// Classify this site as Ts/Tv SNP if biallelic and eligible.
    fn snp_class(&self) -> Option<bool> {
        if !self.is_biallelic {
            return None;
        }
        classify_snp(self.r#ref.as_bytes(), self.first_alt.as_bytes())
    }

    /// Compute total ALT allele count across all sample GT fields.
    fn alt_allele_count(&self) -> u32 {
        let gt_col_idx = self
            .format
            .and_then(|f| f.split(':').position(|t| t == "GT"));
        let Some(gt_idx) = gt_col_idx else {
            return 0;
        };

        // samples_tail holds everything after FORMAT, tab-separated
        self.samples_tail
            .split('\t')
            .map(|sample| {
                let gt = sample.split(':').nth(gt_idx).unwrap_or(".");
                count_alt_alleles(gt)
            })
            .sum()
    }
}

/// Count number of sample columns from the VCF header line.
fn count_samples_from_header(header: &str) -> usize {
    // #CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT [samples...]]
    let cols: Vec<&str> = header.split('\t').collect();
    if cols.len() <= 9 { 0 } else { cols.len() - 9 }
}

/// Compute `--TsTv-by-count` rows.
///
/// Only biallelic SNPs contribute. Max bin = 2 * N_samples - 1 (diploid assumption).
pub fn tstv_by_count<R: BufRead>(reader: R) -> Result<Vec<CountRow>> {
    let mut bins: Vec<(u64, u64)> = Vec::new(); // (n_ts, n_tv) per alt count

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;

        if line.starts_with("##") {
            continue;
        }
        if line.starts_with("#CHROM") {
            let n_samples = count_samples_from_header(&line);
            // vcftools emits bins 0..=2*N_samples-1 (diploid max), total 2*N_samples bins.
            bins = vec![(0, 0); 2 * n_samples];
            continue;
        }

        let Some(rec) = parse_line(&line)? else {
            continue;
        };

        let Some(is_ts) = rec.snp_class() else {
            continue; // indel, multiallelic, or non-SNP — skip entirely
        };

        let count = rec.alt_allele_count() as usize;
        if count < bins.len() {
            if is_ts {
                bins[count].0 += 1;
            } else {
                bins[count].1 += 1;
            }
        }
    }

    let rows = bins
        .into_iter()
        .enumerate()
        .map(|(i, (n_ts, n_tv))| CountRow {
            alt_allele_count: i as u32,
            n_ts,
            n_tv,
            ts_tv: fmt_tstv(n_ts, n_tv),
        })
        .collect();

    Ok(rows)
}

/// Compute `--TsTv-by-qual` rows.
///
/// Thresholds are the unique sorted QUAL values of biallelic SNPs. For each threshold t:
///
/// - LT: sites with QUAL < t
/// - GT: sites with QUAL > t
///
/// QUAL=. is treated as -1.0.
pub fn tstv_by_qual<R: BufRead>(reader: R) -> Result<Vec<QualRow>> {
    // Collect all biallelic SNP sites: (qual, is_ts)
    let mut sites: Vec<(f64, bool)> = Vec::new();

    for line in reader.lines() {
        let line = line.map_err(RsomicsError::Io)?;
        if line.starts_with('#') {
            continue;
        }

        let Some(rec) = parse_line(&line)? else {
            continue;
        };

        let Some(is_ts) = rec.snp_class() else {
            continue;
        };

        // QUAL is always Some here (we parsed it); "." becomes -1.0
        if let Some(q) = rec.qual {
            sites.push((q, is_ts));
        }
    }

    // Collect unique QUAL thresholds from SNP sites, sorted ascending by f64 total order.
    // f64::to_bits() gives wrong ordering for negative values, so we sort by total_cmp.
    let unique_quals: Vec<f64> = {
        let mut vals: Vec<f64> = sites.iter().map(|&(q, _)| q).collect();
        vals.sort_by(f64::total_cmp);
        vals.dedup_by(|a, b| a.to_bits() == b.to_bits());
        vals
    };

    let mut rows = Vec::with_capacity(unique_quals.len());

    for &threshold in &unique_quals {
        let mut n_ts_lt = 0u64;
        let mut n_tv_lt = 0u64;
        let mut n_ts_gt = 0u64;
        let mut n_tv_gt = 0u64;

        for &(q, is_ts) in &sites {
            if q < threshold {
                if is_ts {
                    n_ts_lt += 1;
                } else {
                    n_tv_lt += 1;
                }
            } else if q > threshold {
                if is_ts {
                    n_ts_gt += 1;
                } else {
                    n_tv_gt += 1;
                }
            }
            // q == threshold: excluded from both sides (vcftools semantics)
        }

        let qual_str = format_qual(threshold);
        rows.push(QualRow {
            qual_threshold: qual_str,
            n_ts_lt,
            n_tv_lt,
            ts_tv_lt: fmt_tstv(n_ts_lt, n_tv_lt),
            n_ts_gt,
            n_tv_gt,
            ts_tv_gt: fmt_tstv(n_ts_gt, n_tv_gt),
        });
    }

    Ok(rows)
}

/// Format a QUAL threshold value: integer-looking values print without decimal point.
///
/// vcftools reads QUAL as double and prints with `%g` (which strips trailing zeros).
fn format_qual(q: f64) -> String {
    if q == -1.0 {
        return "-1".to_owned();
    }
    format_g6(q)
}

/// Write `--TsTv-by-count` table to `out`.
pub fn write_count_table<W: std::io::Write>(rows: &[CountRow], mut out: W) -> std::io::Result<()> {
    writeln!(out, "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv")?;
    for r in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            r.alt_allele_count, r.n_ts, r.n_tv, r.ts_tv
        )?;
    }
    Ok(())
}

/// Write `--TsTv-by-qual` table to `out`.
pub fn write_qual_table<W: std::io::Write>(rows: &[QualRow], mut out: W) -> std::io::Result<()> {
    writeln!(
        out,
        "QUAL_THRESHOLD\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD"
    )?;
    for r in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.qual_threshold, r.n_ts_lt, r.n_tv_lt, r.ts_tv_lt, r.n_ts_gt, r.n_tv_gt, r.ts_tv_gt,
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_snp_transitions() {
        assert_eq!(classify_snp(b"A", b"G"), Some(true));
        assert_eq!(classify_snp(b"G", b"A"), Some(true));
        assert_eq!(classify_snp(b"C", b"T"), Some(true));
        assert_eq!(classify_snp(b"T", b"C"), Some(true));
    }

    #[test]
    fn classify_snp_transversions() {
        assert_eq!(classify_snp(b"A", b"C"), Some(false));
        assert_eq!(classify_snp(b"A", b"T"), Some(false));
        assert_eq!(classify_snp(b"G", b"C"), Some(false));
        assert_eq!(classify_snp(b"G", b"T"), Some(false));
        assert_eq!(classify_snp(b"C", b"A"), Some(false));
        assert_eq!(classify_snp(b"T", b"A"), Some(false));
    }

    #[test]
    fn classify_snp_indel_excluded() {
        assert_eq!(classify_snp(b"A", b"AT"), None);
        assert_eq!(classify_snp(b"AT", b"A"), None);
        assert_eq!(classify_snp(b"A", b"<DEL>"), None);
    }

    #[test]
    fn fmt_tstv_cases() {
        assert_eq!(fmt_tstv(0, 0), "nan");
        assert_eq!(fmt_tstv(2, 0), "inf");
        assert_eq!(fmt_tstv(2, 2), "1");
        assert_eq!(fmt_tstv(3, 2), "1.5");
        assert_eq!(fmt_tstv(2, 3), "0.666667");
        assert_eq!(fmt_tstv(7, 3), "2.33333");
    }

    #[test]
    fn format_g6_scientific_and_fixed() {
        assert_eq!(format_g6(1_000_000.0), "1e+06");
        assert_eq!(format_g6(12_345_678.0), "1.23457e+07");
        assert_eq!(format_g6(0.00005), "5e-05");
        assert_eq!(format_g6(0.0001), "0.0001");
        assert_eq!(format_g6(1.0 / 30000.0), "3.33333e-05");
        assert_eq!(format_g6(999999.0), "999999");
        assert_eq!(format_g6(0.0), "0");
    }

    #[test]
    fn count_alt_alleles_cases() {
        assert_eq!(count_alt_alleles("0/1"), 1);
        assert_eq!(count_alt_alleles("1/1"), 2);
        assert_eq!(count_alt_alleles("0/0"), 0);
        assert_eq!(count_alt_alleles("./."), 0);
        assert_eq!(count_alt_alleles("."), 0);
        assert_eq!(count_alt_alleles("0|1"), 1);
        assert_eq!(count_alt_alleles("1|0"), 1);
        assert_eq!(count_alt_alleles("1|1"), 2);
    }
}
