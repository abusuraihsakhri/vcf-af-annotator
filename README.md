# VCF Allele Frequency Annotator

Parse VCF files, extract/calculate allele frequencies, compute population genetics metrics, and classify variants.

## Features

- **VCF Parsing**: Full VCF format support (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLE)
- **Allele Frequency**: Extract AF from INFO tag, calculate from GT/AD fields
- **Hardy-Weinberg Test**: Chi-square test for HWE deviation with p-value calculation
- **Inbreeding Coefficient**: F = 1 - (Ho/He)
- **Variant Classification**: Common (≥1%), Rare (0.1-1%), Ultra-rare (<0.1%)
- **gnomAD Comparison**: Fold-difference against population frequencies

## Installation

Zero dependencies — Python 3.8+ stdlib only.

```bash
cd vcf-af-annotator
```

## Usage

### CLI

```bash
# Annotate a VCF file
python cli.py annotate -i input.vcf -o annotated.csv

# Hardy-Weinberg test
python cli.py hwe --hom-ref 25 --het 50 --hom-alt 25

# Classify a variant
python cli.py classify 0.05

# Parse a single VCF line
python cli.py parse "chr1	100	rs123	A	G	30	PASS	AF=0.25"

# Batch process from CSV
python cli.py batch -i variants.csv -o results.csv
```

### Python API

```python
from vcf_af import (
    parse_vcf_line,
    annotate_vcf_record,
    chi_square_hwe,
    classify_variant_by_af,
    hardy_weinberg_expected,
)

# Parse a VCF line
rec = parse_vcf_line("chr1\t100\trs123\tA\tG\t30\tPASS\tAF=0.25")

# Annotate with AF and HWE
ann = annotate_vcf_record(rec)
print(f"AF: {ann['AF']}, Category: {ann['category']}")

# Hardy-Weinberg test
result = chi_square_hwe(25, 50, 25, 0.5)
print(f"Chi²: {result['chi_square']}, p={result['p_value']}, In HWE: {result['in_hwe']}")

# Expected genotype frequencies
hw = hardy_weinberg_expected(0.3)
print(f"p²={hw['p_squared']:.3f}, 2pq={hw['two_pq']:.3f}, q²={hw['q_squared']:.3f}")
```

## Formulas

### Hardy-Weinberg Expected Frequencies
- p² (homozygous reference)
- 2pq (heterozygous)
- q² (homozygous alternate)

### Chi-Square Test
χ² = Σ((observed - expected)² / expected), df = 1

### Inbreeding Coefficient
F = 1 - (Ho / He)

## Tests

```bash
python -m pytest test_vcf_af.py -v
```

## License

MIT
