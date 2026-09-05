#!/usr/bin/env python3
"""
VCF Allele Frequency Annotator — Real Implementation

Parses VCF files, extracts/calculates allele frequencies, computes population
genetics metrics (Hardy-Weinberg, chi-square), and classifies variants by AF.

Zero external dependencies — Python stdlib only.
"""

import csv
import math
import os
import re
import sys
from typing import Dict, Any, List, Optional, Tuple


# ---------------------------------------------------------------------------
# VCF Parsing
# ---------------------------------------------------------------------------

def parse_vcf_line(line: str) -> Optional[Dict[str, Any]]:
    """Parse a single VCF data line into a dict.

    VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT SAMPLE...]
    """
    line = line.strip()
    if not line or line.startswith('#'):
        return None

    fields = line.split('\t')
    if len(fields) < 8:
        return None

    record = {
        'CHROM': fields[0],
        'POS': int(fields[1]),
        'ID': fields[2],
        'REF': fields[3],
        'ALT': fields[4],
        'QUAL': fields[5] if fields[5] != '.' else None,
        'FILTER': fields[6],
        'INFO': fields[7],
    }

    # Parse INFO field into dict
    record['INFO_DICT'] = parse_info(fields[7])

    # Parse genotype fields if present
    if len(fields) >= 10:
        record['FORMAT'] = fields[8]
        record['SAMPLE'] = fields[9]
        record['GENOTYPE'] = parse_genotype(fields[8], fields[9])

    return record


def parse_info(info_str: str) -> Dict[str, Any]:
    """Parse VCF INFO field into a dictionary.

    Format: KEY1=val1;KEY2=val2;FLAG
    """
    info = {}
    if not info_str or info_str == '.':
        return info

    for entry in info_str.split(';'):
        if '=' in entry:
            key, val = entry.split('=', 1)
            # Try to parse numeric values
            try:
                if ',' in val:
                    info[key] = [float(v) if '.' in v else int(v) for v in val.split(',')]
                elif '.' in val:
                    info[key] = float(val)
                else:
                    info[key] = int(val)
            except ValueError:
                info[key] = val
        else:
            info[entry] = True  # Flag

    return info


def parse_genotype(format_str: str, sample_str: str) -> Dict[str, Any]:
    """Parse VCF FORMAT/SAMPLE fields into a dictionary.

    Handles GT, DP, AD, GQ, PL fields.
    """
    fmt_keys = format_str.split(':')
    sample_vals = sample_str.split(':')
    gt = {}

    for i, key in enumerate(fmt_keys):
        if i < len(sample_vals):
            val = sample_vals[i]
            if key == 'GT':
                gt['GT'] = val
                gt['alleles'] = parse_gt_alleles(val)
            elif key == 'DP':
                try:
                    gt['DP'] = int(val)
                except ValueError:
                    gt['DP'] = val
            elif key == 'AD':
                try:
                    gt['AD'] = [int(v) for v in val.split(',')]
                except ValueError:
                    gt['AD'] = val
            elif key == 'GQ':
                try:
                    gt['GQ'] = int(val)
                except ValueError:
                    gt['GQ'] = val
            elif key == 'PL':
                try:
                    gt['PL'] = [int(v) for v in val.split(',')]
                except ValueError:
                    gt['PL'] = val
            else:
                gt[key] = val

    return gt


def parse_gt_alleles(gt_str: str) -> List[int]:
    """Parse genotype string like '0/1' or '1|0' into allele indices."""
    # Handle both / and | separators
    sep = '/' if '/' in gt_str else '|'
    parts = gt_str.split(sep)
    alleles = []
    for p in parts:
        try:
            alleles.append(int(p))
        except ValueError:
            alleles.append(-1)  # Missing
    return alleles


# ---------------------------------------------------------------------------
# Allele Frequency Extraction
# ---------------------------------------------------------------------------

def extract_af_from_info(info_dict: Dict[str, Any]) -> Optional[float]:
    """Extract allele frequency from VCF INFO field (AF tag)."""
    af = info_dict.get('AF')
    if af is not None:
        if isinstance(af, list):
            return af[0]  # First ALT allele frequency
        return float(af)
    return None


def calculate_af_from_genotypes(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Calculate allele frequency from genotype (GT) fields across samples.

    Returns dict with allele counts, frequencies, and genotype counts.
    """
    alt_count = 0
    total_alleles = 0
    genotype_counts = {'hom_ref': 0, 'het': 0, 'hom_alt': 0, 'missing': 0}

    for rec in records:
        gt = rec.get('GENOTYPE', {})
        alleles = gt.get('alleles', [])
        if not alleles:
            genotype_counts['missing'] += 1
            continue

        for a in alleles:
            if a < 0:
                genotype_counts['missing'] += 1
                continue
            total_alleles += 1
            if a > 0:
                alt_count += 1

        if len(alleles) == 2:
            if alleles[0] == 0 and alleles[1] == 0:
                genotype_counts['hom_ref'] += 1
            elif alleles[0] > 0 and alleles[1] > 0:
                genotype_counts['hom_alt'] += 1
            elif alleles[0] != alleles[1]:
                genotype_counts['het'] += 1

    af = alt_count / total_alleles if total_alleles > 0 else 0.0

    return {
        'alt_count': alt_count,
        'total_alleles': total_alleles,
        'af': af,
        'ref_freq': 1.0 - af,
        'genotype_counts': genotype_counts,
    }


def calculate_af_from_ad(info_dict: Dict[str, Any],
                          genotype: Dict[str, Any]) -> Optional[float]:
    """Calculate allele frequency from AD (allelic depth) field."""
    ad = genotype.get('AD')
    if ad and isinstance(ad, list) and len(ad) >= 2:
        total = sum(ad)
        if total > 0:
            return ad[1] / total  # ALT depth / total depth
    return None


# ---------------------------------------------------------------------------
# Population Genetics Metrics
# ---------------------------------------------------------------------------

def hardy_weinberg_expected(af: float) -> Dict[str, float]:
    """Calculate expected Hardy-Weinberg genotype frequencies.

    Given allele frequency p (ref) and q (alt):
        p² (homozygous reference)
        2pq (heterozygous)
        q² (homozygous alternate)
    """
    p = 1.0 - af  # Reference allele frequency
    q = af         # Alternate allele frequency

    return {
        'p': p,
        'q': q,
        'p_squared': p * p,      # Expected hom-ref frequency
        'two_pq': 2.0 * p * q,   # Expected het frequency
        'q_squared': q * q,       # Expected hom-alt frequency
    }


def chi_square_hwe(observed_hom_ref: int, observed_het: int,
                    observed_hom_alt: int, af: Optional[float] = None) -> Dict[str, Any]:
    """Chi-square test for Hardy-Weinberg Equilibrium deviation.

    If af is not provided, it is calculated from observed genotype counts.

    Returns:
        Dict with chi_square statistic, p_value, degrees_of_freedom, in_hwe.
    """
    n = observed_hom_ref + observed_het + observed_hom_alt
    if n == 0:
        return {'chi_square': 0.0, 'p_value': 1.0, 'df': 1, 'in_hwe': True,
                'error': 'No genotypes observed'}

    # Calculate allele frequency from counts if not provided
    if af is None:
        total_alleles = 2 * n
        alt_alleles = 2 * observed_hom_alt + observed_het
        af = alt_alleles / total_alleles

    hw = hardy_weinberg_expected(af)

    # Expected counts
    exp_hom_ref = hw['p_squared'] * n
    exp_het = hw['two_pq'] * n
    exp_hom_alt = hw['q_squared'] * n

    # Chi-square calculation
    chi_sq = 0.0
    for obs, exp in [(observed_hom_ref, exp_hom_ref),
                     (observed_het, exp_het),
                     (observed_hom_alt, exp_hom_alt)]:
        if exp > 0:
            chi_sq += (obs - exp) ** 2 / exp

    # P-value approximation using incomplete gamma function approximation
    df = 1  # One degree of freedom for HWE test
    p_value = chi_square_p_value(chi_sq, df)

    return {
        'chi_square': round(chi_sq, 4),
        'p_value': p_value,
        'df': df,
        'in_hwe': p_value > 0.05,
        'expected_hom_ref': round(exp_hom_ref, 2),
        'expected_het': round(exp_het, 2),
        'expected_hom_alt': round(exp_hom_alt, 2),
        'observed_hom_ref': observed_hom_ref,
        'observed_het': observed_het,
        'observed_hom_alt': observed_hom_alt,
        'sample_size': n,
        'allele_freq_alt': round(af, 6),
    }


def chi_square_p_value(chi_sq: float, df: int) -> float:
    """Approximate p-value for chi-square distribution using regularized
    incomplete gamma function (stdlib only).

    For df=1: P(X > x) = 2 * (1 - Phi(sqrt(x))) where Phi is normal CDF.
    We use math.erfc for the normal CDF approximation.
    """
    if chi_sq <= 0:
        return 1.0
    if df == 1:
        # For 1 df: p = 2 * (1 - Phi(sqrt(chi_sq)))
        # Phi(x) = 0.5 * (1 + erf(x / sqrt(2)))
        # So p = erfc(sqrt(chi_sq / 2))
        return math.erfc(math.sqrt(chi_sq / 2.0))
    elif df == 2:
        # For 2 df: p = exp(-chi_sq / 2)
        return math.exp(-chi_sq / 2.0)
    else:
        # General case: use Wilson-Hilferty approximation
        # For large df, (chi_sq/df)^(1/3) ~ Normal(1 - 2/(9*df), 2/(9*df))
        k = df
        x = chi_sq
        # Regularized incomplete gamma: P(k/2, x/2) = 1 - p_value
        # Use series expansion for small x, continued fraction for large x
        a = k / 2.0
        z = x / 2.0
        if z < a + 1:
            return 1.0 - _gammainc_series(a, z)
        else:
            return _gammainc_cf(a, z)


def _gammainc_series(a: float, x: float) -> float:
    """Series expansion for regularized lower incomplete gamma function."""
    if x == 0:
        return 0.0
    ap = a
    total = 1.0 / a
    term = 1.0 / a
    for _ in range(300):
        ap += 1.0
        term *= x / ap
        total += term
        if abs(term) < abs(total) * 1e-15:
            break
    return total * math.exp(-x + a * math.log(x) - math.lgamma(a))


def _gammainc_cf(a: float, x: float) -> float:
    """Continued fraction for regularized upper incomplete gamma function."""
    b = x + 1.0 - a
    c = 1e30
    d = 1.0 / b
    h = d
    for i in range(1, 300):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < 1e-30:
            d = 1e-30
        c = b + an / c
        if abs(c) < 1e-30:
            c = 1e-30
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-15:
            break
    return h * math.exp(-x + a * math.log(x) - math.lgamma(a))


# ---------------------------------------------------------------------------
# Inbreeding Coefficient
# ---------------------------------------------------------------------------

def inbreeding_coefficient(observed_het: int, expected_het: float) -> Dict[str, Any]:
    """Calculate inbreeding coefficient F = 1 - (Ho / He).

    F > 0: excess homozygosity (inbreeding)
    F < 0: excess heterozygosity
    F = 0: HWE
    """
    if expected_het == 0:
        return {'F': 0.0, 'interpretation': 'No heterozygosity expected',
                'observed_het': observed_het, 'expected_het': expected_het}

    ho = observed_het
    he = expected_het
    f = 1.0 - (ho / he)

    if f > 0.1:
        interp = 'Excess homozygosity — possible inbreeding or population structure'
    elif f < -0.1:
        interp = 'Excess heterozygosity — possible outbreeding or selection against homozygotes'
    else:
        interp = 'Consistent with Hardy-Weinberg Equilibrium'

    return {
        'F': round(f, 4),
        'observed_heterozygosity': ho,
        'expected_heterozygosity': round(he, 4),
        'interpretation': interp,
    }


# ---------------------------------------------------------------------------
# Variant Classification
# ---------------------------------------------------------------------------

def classify_variant_by_af(af: float) -> Dict[str, Any]:
    """Classify a variant by allele frequency.

    Categories:
        Common:     AF >= 1%     (polymorphism)
        Rare:       0.1% <= AF < 1%
        Ultra-rare: AF < 0.1%
    """
    if af >= 0.01:
        category = 'common'
        description = 'Common variant / polymorphism (AF >= 1%)'
    elif af >= 0.001:
        category = 'rare'
        description = 'Rare variant (0.1% <= AF < 1%)'
    else:
        category = 'ultra-rare'
        description = 'Ultra-rare variant (AF < 0.1%)'

    return {
        'af': af,
        'category': category,
        'description': description,
    }


def compare_to_gnomad(af: float, gnomad_af: Optional[float] = None) -> Dict[str, Any]:
    """Compare observed AF to gnomAD-style population frequency.

    If gnomad_af is provided, computes fold-difference.
    Otherwise, classifies against standard gnomAD-like thresholds.
    """
    result = classify_variant_by_af(af)

    if gnomad_af is not None:
        if gnomad_af > 0:
            fold_diff = af / gnomad_af
            result['gnomad_af'] = gnomad_af
            result['fold_difference'] = round(fold_diff, 4)
            if fold_diff > 2.0:
                result['comparison'] = 'Observed AF significantly higher than gnomAD'
            elif fold_diff < 0.5:
                result['comparison'] = 'Observed AF significantly lower than gnomAD'
            else:
                result['comparison'] = 'Observed AF consistent with gnomAD'
        else:
            result['gnomad_af'] = 0.0
            result['comparison'] = 'Variant absent from gnomAD'

    return result


# ---------------------------------------------------------------------------
# Full VCF Annotation
# ---------------------------------------------------------------------------

def annotate_vcf_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Annotate a single VCF record with AF and population genetics metrics."""
    info = record.get('INFO_DICT', {})
    genotype = record.get('GENOTYPE', {})

    # Try to get AF from INFO first
    af = extract_af_from_info(info)

    # Fall back to AD field
    if af is None and genotype:
        af = calculate_af_from_ad(info, genotype)

    # Fall back to GT field
    if af is None and genotype:
        alleles = genotype.get('alleles', [])
        if alleles:
            alt_count = sum(1 for a in alleles if a > 0)
            af = alt_count / len(alleles) if alleles else 0.0

    if af is None:
        af = 0.0

    # Classify
    classification = classify_variant_by_af(af)

    # HWE (from genotype if available)
    hwe_result = None
    if genotype and genotype.get('alleles'):
        alleles = genotype['alleles']
        if len(alleles) == 2:
            hom_ref = 1 if alleles == [0, 0] else 0
            het = 1 if len(set(alleles)) > 1 and all(a >= 0 for a in alleles) else 0
            hom_alt = 1 if alleles == [1, 1] else 0
            hwe_result = chi_square_hwe(hom_ref, het, hom_alt, af)

    result = {
        'CHROM': record.get('CHROM'),
        'POS': record.get('POS'),
        'REF': record.get('REF'),
        'ALT': record.get('ALT'),
        'AF': round(af, 6),
        'category': classification['category'],
        'description': classification['description'],
    }

    if hwe_result:
        result['hwe_chi_square'] = hwe_result['chi_square']
        result['hwe_p_value'] = hwe_result['p_value']
        result['in_hwe'] = hwe_result['in_hwe']

    return result


def _validate_filepath(filepath: str) -> str:
    """Validate filepath for path traversal and null byte injection."""
    if not filepath:
        raise ValueError("Filepath must not be empty")
    if '\x00' in filepath:
        raise ValueError("Filepath contains null bytes")
    # Resolve to absolute path and check for traversal attempts
    normalized = os.path.normpath(filepath)
    if '..' in normalized.split(os.sep):
        raise ValueError(f"Path traversal detected in: {filepath}")
    return normalized


def parse_vcf_file(filepath: str) -> List[Dict[str, Any]]:
    """Parse a VCF file and return list of record dicts."""
    filepath = _validate_filepath(filepath)
    records = []
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue
            rec = parse_vcf_line(line)
            if rec:
                records.append(rec)
    return records


def annotate_vcf_file(input_path: str, output_path: str) -> List[Dict[str, Any]]:
    """Annotate an entire VCF file and write results to CSV."""
    input_path = _validate_filepath(input_path)
    output_path = _validate_filepath(output_path)
    records = parse_vcf_file(input_path)
    annotations = [annotate_vcf_record(r) for r in records]

    if annotations:
        fieldnames = list(annotations[0].keys())
        with open(output_path, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(annotations)

    return annotations


def process_batch(input_csv: str, output_csv: str) -> List[Dict[str, Any]]:
    """Process a CSV of variant data with AF values.

    Expected columns: CHROM, POS, REF, ALT, AF (or INFO field).
    """
    input_csv = _validate_filepath(input_csv)
    output_csv = _validate_filepath(output_csv)
    with open(input_csv, mode='r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)

    out_fields = fieldnames + ['af_category', 'hwe_p_value', 'in_hwe']
    out_rows = []

    for r in rows:
        af = None
        # Try AF column
        for key in ['AF', 'af', 'allele_frequency']:
            if key in r and r[key]:
                try:
                    af = float(r[key])
                    break
                except ValueError:
                    pass

        # Try INFO column
        if af is None and 'INFO' in r:
            info = parse_info(r['INFO'])
            af = extract_af_from_info(info)

        if af is None:
            af = 0.0

        classification = classify_variant_by_af(af)
        row_dict = dict(r)
        row_dict['af_category'] = classification['category']
        out_rows.append(row_dict)

    with open(output_csv, mode='w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=out_fields)
        writer.writeheader()
        writer.writerows(out_rows)

    return out_rows
