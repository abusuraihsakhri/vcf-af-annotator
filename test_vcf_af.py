#!/usr/bin/env python3
"""Tests for VCF Allele Frequency Annotator — 20+ real functional tests."""
import os

import pytest
from vcf_af import (
    parse_vcf_line,
    parse_info,
    parse_genotype,
    parse_gt_alleles,
    extract_af_from_info,
    calculate_af_from_genotypes,
    calculate_af_from_ad,
    hardy_weinberg_expected,
    chi_square_hwe,
    chi_square_p_value,
    inbreeding_coefficient,
    classify_variant_by_af,
    compare_to_gnomad,
    annotate_vcf_record,
    process_batch,
    _validate_filepath,
)


# ---------------------------------------------------------------------------
# VCF Parsing Tests
# ---------------------------------------------------------------------------

class TestVCFParsing:
    def test_parse_simple_line(self):
        line = "chr1\t12345\trs123\tA\tG\t30\tPASS\tAF=0.25"
        rec = parse_vcf_line(line)
        assert rec is not None
        assert rec['CHROM'] == 'chr1'
        assert rec['POS'] == 12345
        assert rec['REF'] == 'A'
        assert rec['ALT'] == 'G'

    def test_parse_header_returns_none(self):
        assert parse_vcf_line('#CHROM\tPOS\tID') is None

    def test_parse_empty_returns_none(self):
        assert parse_vcf_line('') is None

    def test_parse_with_genotype(self):
        line = "chr1\t100\t.\tA\tG\t50\tPASS\tAF=0.5\tGT:DP\t0/1:30"
        rec = parse_vcf_line(line)
        assert rec['GENOTYPE']['GT'] == '0/1'
        assert rec['GENOTYPE']['DP'] == 30

    def test_parse_info_field(self):
        info = parse_info('AF=0.25;DP=100;DB')
        assert info['AF'] == 0.25
        assert info['DP'] == 100
        assert info['DB'] is True

    def test_parse_info_comma_separated(self):
        info = parse_info('AF=0.1,0.05')
        assert isinstance(info['AF'], list)
        assert len(info['AF']) == 2

    def test_parse_info_dot(self):
        info = parse_info('.')
        assert len(info) == 0

    def test_parse_gt_slash(self):
        alleles = parse_gt_alleles('0/1')
        assert alleles == [0, 1]

    def test_parse_gt_pipe(self):
        alleles = parse_gt_alleles('1|0')
        assert alleles == [1, 0]

    def test_parse_gt_homozygous(self):
        alleles = parse_gt_alleles('1/1')
        assert alleles == [1, 1]


# ---------------------------------------------------------------------------
# Allele Frequency Tests
# ---------------------------------------------------------------------------

class TestAlleleFrequency:
    def test_extract_af_from_info(self):
        info = {'AF': 0.25, 'DP': 100}
        assert extract_af_from_info(info) == 0.25

    def test_extract_af_list(self):
        info = {'AF': [0.25, 0.1]}
        assert extract_af_from_info(info) == 0.25

    def test_extract_af_missing(self):
        info = {'DP': 100}
        assert extract_af_from_info(info) is None

    def test_calculate_af_from_genotypes(self):
        records = [
            {'GENOTYPE': {'alleles': [0, 0]}},
            {'GENOTYPE': {'alleles': [0, 1]}},
            {'GENOTYPE': {'alleles': [1, 1]}},
        ]
        result = calculate_af_from_genotypes(records)
        # 3 samples, 6 alleles, 3 alt alleles -> AF = 0.5
        assert result['af'] == 0.5
        assert result['alt_count'] == 3
        assert result['total_alleles'] == 6

    def test_calculate_af_all_ref(self):
        records = [
            {'GENOTYPE': {'alleles': [0, 0]}},
            {'GENOTYPE': {'alleles': [0, 0]}},
        ]
        result = calculate_af_from_genotypes(records)
        assert result['af'] == 0.0

    def test_calculate_af_from_ad(self):
        info = {}
        genotype = {'AD': [10, 30]}
        af = calculate_af_from_ad(info, genotype)
        assert af == 0.75  # 30 / (10+30)


# ---------------------------------------------------------------------------
# Hardy-Weinberg Tests
# ---------------------------------------------------------------------------

class TestHardyWeinberg:
    def test_hw_expected_frequencies(self):
        hw = hardy_weinberg_expected(0.5)
        assert hw['p'] == 0.5
        assert hw['q'] == 0.5
        assert abs(hw['p_squared'] - 0.25) < 1e-10
        assert abs(hw['two_pq'] - 0.5) < 1e-10
        assert abs(hw['q_squared'] - 0.25) < 1e-10

    def test_hw_frequencies_sum_to_one(self):
        hw = hardy_weinberg_expected(0.3)
        total = hw['p_squared'] + hw['two_pq'] + hw['q_squared']
        assert abs(total - 1.0) < 1e-10

    def test_hw_rare_allele(self):
        hw = hardy_weinberg_expected(0.01)
        assert hw['q_squared'] < 0.001  # Very rare homozygous alt

    def test_chi_square_hwe_perfect(self):
        # At AF=0.5, expected is 25:50:25 for n=100
        result = chi_square_hwe(25, 50, 25, 0.5)
        assert result['chi_square'] < 0.01
        assert result['in_hwe'] is True

    def test_chi_square_hwe_deviation(self):
        # Excess homozygosity
        result = chi_square_hwe(40, 20, 40, 0.5)
        assert result['chi_square'] > 10.0
        assert result['in_hwe'] is False

    def test_chi_square_auto_af(self):
        # AF calculated from counts: (2*25 + 50) / 200 = 0.5
        result = chi_square_hwe(25, 50, 25)
        assert result['in_hwe'] is True

    def test_chi_square_zero_samples(self):
        result = chi_square_hwe(0, 0, 0)
        assert 'error' in result


# ---------------------------------------------------------------------------
# Chi-Square P-Value Tests
# ---------------------------------------------------------------------------

class TestChiSquarePValue:
    def test_zero_chi_square(self):
        assert chi_square_p_value(0.0, 1) == 1.0

    def test_large_chi_square(self):
        p = chi_square_p_value(10.0, 1)
        assert p < 0.01  # Very significant

    def test_small_chi_square(self):
        p = chi_square_p_value(0.5, 1)
        assert p > 0.05  # Not significant

    def test_two_degrees_freedom(self):
        p = chi_square_p_value(6.0, 2)
        assert 0.01 < p < 0.1


# ---------------------------------------------------------------------------
# Inbreeding Coefficient Tests
# ---------------------------------------------------------------------------

class TestInbreedingCoefficient:
    def test_no_inbreeding(self):
        result = inbreeding_coefficient(50, 50.0)
        assert result['F'] == 0.0

    def test_inbreeding(self):
        result = inbreeding_coefficient(20, 50.0)
        assert result['F'] == 0.6  # 1 - 20/50
        assert 'inbreeding' in result['interpretation'].lower() or 'homozygosity' in result['interpretation'].lower()

    def test_excess_het(self):
        result = inbreeding_coefficient(60, 50.0)
        assert result['F'] < 0


# ---------------------------------------------------------------------------
# Variant Classification Tests
# ---------------------------------------------------------------------------

class TestVariantClassification:
    def test_common_variant(self):
        result = classify_variant_by_af(0.05)
        assert result['category'] == 'common'

    def test_rare_variant(self):
        result = classify_variant_by_af(0.005)
        assert result['category'] == 'rare'

    def test_ultra_rare_variant(self):
        result = classify_variant_by_af(0.0005)
        assert result['category'] == 'ultra-rare'

    def test_boundary_common_rare(self):
        result = classify_variant_by_af(0.01)
        assert result['category'] == 'common'

    def test_boundary_rare_ultra(self):
        result = classify_variant_by_af(0.001)
        assert result['category'] == 'rare'


# ---------------------------------------------------------------------------
# gnomAD Comparison Tests
# ---------------------------------------------------------------------------

class TestGnomADComparison:
    def test_consistent(self):
        result = compare_to_gnomad(0.05, 0.04)
        assert 'consistent' in result['comparison'].lower()

    def test_higher_than_gnomad(self):
        result = compare_to_gnomad(0.10, 0.01)
        assert 'higher' in result['comparison'].lower()

    def test_without_gnomad(self):
        result = compare_to_gnomad(0.05)
        assert 'category' in result


# ---------------------------------------------------------------------------
# Full Annotation Tests
# ---------------------------------------------------------------------------

class TestFullAnnotation:
    def test_annotate_record(self):
        line = "chr1\t100\trs123\tA\tG\t30\tPASS\tAF=0.25\tGT:DP\t0/1:30"
        rec = parse_vcf_line(line)
        ann = annotate_vcf_record(rec)
        assert ann['AF'] == 0.25
        assert ann['category'] == 'common'

    def test_batch_csv(self, tmp_path):
        csv_in = tmp_path / 'in.csv'
        csv_out = tmp_path / 'out.csv'
        csv_in.write_text(
            'CHROM,POS,REF,ALT,AF\n'
            'chr1,100,A,G,0.05\n'
            'chr2,200,C,T,0.0005\n',
            encoding='utf-8'
        )
        results = process_batch(str(csv_in), str(csv_out))
        assert len(results) == 2
        assert csv_out.exists()


# ---------------------------------------------------------------------------
# Filepath Validation Tests
# ---------------------------------------------------------------------------

class TestFilepathValidation:
    def test_valid_filepath(self):
        result = _validate_filepath("data/sample.vcf")
        assert result == os.path.normpath("data/sample.vcf")

    def test_empty_filepath_raises(self):
        with pytest.raises(ValueError, match="empty"):
            _validate_filepath("")

    def test_null_byte_raises(self):
        with pytest.raises(ValueError, match="null bytes"):
            _validate_filepath("file\x00name.vcf")

    def test_path_traversal_raises(self):
        with pytest.raises(ValueError, match="traversal"):
            _validate_filepath("../../../etc/passwd")

    def test_path_traversal_with_filename_raises(self):
        with pytest.raises(ValueError, match="traversal"):
            _validate_filepath("data/../../../etc/shadow")


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
