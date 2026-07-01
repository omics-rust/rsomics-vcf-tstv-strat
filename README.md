# rsomics-vcf-tstv-strat

Stratified transition/transversion (Ts/Tv) statistics from a VCF — reimplements
`vcftools --TsTv-by-count` and `--TsTv-by-qual`.

## Usage

```
rsomics-vcf-tstv-strat [OPTIONS] [INPUT]
```

```
OPTIONS:
  -b, --by <MODE>   count  — stratify by ALT allele count across samples (default)
                    qual   — stratify by QUAL threshold
```

### By allele count

```sh
rsomics-vcf-tstv-strat --by count input.vcf
rsomics-vcf-tstv-strat input.vcf          # count is the default
bcftools view cohort.bcf | rsomics-vcf-tstv-strat
```

Output columns (tab-separated, one row per allele count from 0 to 2×N_samples−1):

```
ALT_ALLELE_COUNT  N_Ts  N_Tv  Ts/Tv
```

### By QUAL threshold

```sh
rsomics-vcf-tstv-strat --by qual input.vcf
```

Output columns (one row per unique QUAL value observed at biallelic SNP sites):

```
QUAL_THRESHOLD  N_Ts_LT_QUAL_THRESHOLD  N_Tv_LT_QUAL_THRESHOLD  Ts/Tv_LT_QUAL_THRESHOLD
N_Ts_GT_QUAL_THRESHOLD  N_Tv_GT_QUAL_THRESHOLD  Ts/Tv_GT_QUAL_THRESHOLD
```

`LT` counts sites with QUAL strictly less than the threshold; `GT` counts sites
with QUAL strictly greater. Sites at exactly the threshold appear in neither column.

## Semantics

**by-count**

- Only biallelic SNPs (single-base REF, single-base first ALT) are included;
  indels, multiallelic sites, and symbolic alleles are fully excluded.
- `ALT_ALLELE_COUNT` = sum of non-ref, non-missing allele calls in the GT field
  across all samples (`0/1` → 1, `1/1` → 2, `./.` or `.` → 0).
- Bins span `0 .. 2×N_samples` (diploid assumption, fixed regardless of observed
  counts).
- Ts/Tv formatted as C `%g` with precision 6; `nan` when Ts=Tv=0, `inf` when
  Tv=0 and Ts>0.

**by-qual**

- Thresholds are the unique sorted QUAL values of biallelic SNP sites only
  (indel QUALs do not generate rows).
- QUAL `.` (missing) is treated as −1 and creates a row with
  `QUAL_THRESHOLD=-1`.
- Same Ts/Tv classification as by-count (biallelic SNPs only).
- Note: vcftools 0.1.17 has an uninitialized-memory bug in the `N_Tv_GT` column
  of the last threshold row; this implementation outputs the correct value (0)
  instead.

## Origin

This crate is an independent Rust reimplementation based on:

- Black-box behavioral testing against `vcftools 0.1.17`
- The VCF format specification (v4.2)
- The vcftools paper (Danecek et al., 2011, DOI: 10.1093/bioinformatics/btr330)

No GPL source code was used as reference. Fixtures were generated independently.

License: MIT OR Apache-2.0  
Upstream credit: vcftools 0.1.17 <https://github.com/vcftools/vcftools> (GPL-3.0)
