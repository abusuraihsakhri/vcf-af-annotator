#!/usr/bin/env python3
"""
CLI for VCF Allele Frequency Annotator.

Subcommands:
  annotate  — Annotate a VCF file with AF and population genetics
  hwe       — Run Hardy-Weinberg test on genotype counts
  classify  — Classify a variant by allele frequency
  parse     — Parse and display a VCF record
  batch     — Batch-annotate from CSV
"""
import argparse
import json
import sys

from vcf_af import (
    parse_vcf_line,
    annotate_vcf_record,
    annotate_vcf_file,
    chi_square_hwe,
    classify_variant_by_af,
    hardy_weinberg_expected,
    inbreeding_coefficient,
    process_batch,
)


def cmd_annotate(args):
    """Annotate a VCF file."""
    results = annotate_vcf_file(args.input, args.output)
    print(f"Annotated {len(results)} variants -> {args.output}")
    # Print summary
    categories = {}
    for r in results:
        cat = r.get('category', 'unknown')
        categories[cat] = categories.get(cat, 0) + 1
    print(f"Categories: {json.dumps(categories)}")


def cmd_hwe(args):
    """Run Hardy-Weinberg equilibrium test."""
    result = chi_square_hwe(args.hom_ref, args.het, args.hom_alt, args.af)
    print(json.dumps(result, indent=2))


def cmd_classify(args):
    """Classify a variant by allele frequency."""
    result = classify_variant_by_af(args.af)
    print(json.dumps(result, indent=2))


def cmd_parse(args):
    """Parse a single VCF line."""
    record = parse_vcf_line(args.vcf_line)
    if record:
        annotation = annotate_vcf_record(record)
        print(json.dumps(annotation, indent=2, default=str))
    else:
        print("Could not parse VCF line", file=sys.stderr)
        sys.exit(1)


def cmd_batch(args):
    """Batch process from CSV."""
    results = process_batch(args.input, args.output)
    print(f"Processed {len(results)} variants -> {args.output}")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog='vcf-af-annotator',
        description='VCF Allele Frequency Annotator — Parse VCF, compute AF, HWE, classify variants'
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # annotate
    p_ann = subparsers.add_parser('annotate', help='Annotate a VCF file')
    p_ann.add_argument('-i', '--input', required=True, help='Input VCF file')
    p_ann.add_argument('-o', '--output', default='annotated.csv', help='Output CSV file')
    p_ann.set_defaults(func=cmd_annotate)

    # hwe
    p_hwe = subparsers.add_parser('hwe', help='Hardy-Weinberg equilibrium test')
    p_hwe.add_argument('--hom-ref', type=int, required=True, help='Homozygous reference count')
    p_hwe.add_argument('--het', type=int, required=True, help='Heterozygous count')
    p_hwe.add_argument('--hom-alt', type=int, required=True, help='Homozygous alternate count')
    p_hwe.add_argument('--af', type=float, default=None, help='Alternate allele frequency (auto-calculated if omitted)')
    p_hwe.set_defaults(func=cmd_hwe)

    # classify
    p_cls = subparsers.add_parser('classify', help='Classify variant by AF')
    p_cls.add_argument('af', type=float, help='Allele frequency (0-1)')
    p_cls.set_defaults(func=cmd_classify)

    # parse
    p_parse = subparsers.add_parser('parse', help='Parse a single VCF line')
    p_parse.add_argument('vcf_line', help='VCF data line (tab-separated)')
    p_parse.set_defaults(func=cmd_parse)

    # batch
    p_batch = subparsers.add_parser('batch', help='Batch annotate from CSV')
    p_batch.add_argument('-i', '--input', required=True, help='Input CSV file')
    p_batch.add_argument('-o', '--output', default='results.csv', help='Output CSV file')
    p_batch.set_defaults(func=cmd_batch)

    args = parser.parse_args(argv)
    args.func(args)


if __name__ == '__main__':
    main()
