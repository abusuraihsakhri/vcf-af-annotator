#!/usr/bin/env python3
"""
CLI for VCF Allele Frequency Annotator.

Subcommands:
  annotate  — Annotate a VCF file with AF and population genetics
  hwe       — Run Hardy-Weinberg test on genotype counts
  classify  — Classify a variant by allele frequency
  parse     — Parse and display a VCF record
  batch     — Batch-annotate from CSV
  audit     — Dispatch an enterprise audit task via the supervisor
  chat      — Query the air-gapped supervisory assistant
  verify-audit — Verify the HMAC-SHA256 audit trail integrity
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


def cmd_audit(args):
    """Dispatch an enterprise audit task via the supervisor."""
    from agents.supervisor import SystemSupervisor
    from agents.models import SystemTaskPayload
    supervisor = SystemSupervisor(model_provider="mock")
    payload = SystemTaskPayload(
        task_id=args.task_id,
        target_identifier=args.target_identifier,
        primary_metric=args.primary_metric,
        secondary_metric=args.secondary_metric,
        status_descriptor=args.status_descriptor,
        is_critical_flag=args.is_critical_flag,
    )
    dossier = supervisor.process_task(payload)
    print(json.dumps(dossier.to_dict(), indent=2, default=str))


def cmd_chat(args):
    """Query the air-gapped supervisory assistant."""
    from agents.supervisor import SystemSupervisor
    supervisor = SystemSupervisor(model_provider="mock")
    query = " ".join(args.query)
    response = supervisor.query_supervisory_chat(query)
    print(response)


def cmd_verify_audit(args):
    """Verify the HMAC-SHA256 audit trail integrity."""
    from agents.base import AuditLogger
    verified = AuditLogger.verify_integrity()
    trail_len = len(AuditLogger.get_trail())
    status = "VALID" if verified else "TAMPERED"
    print(f"Audit trail: {trail_len} blocks, integrity: {status}")
    if not verified:
        sys.exit(1)


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

    # audit
    p_audit = subparsers.add_parser('audit', help='Dispatch enterprise audit task via supervisor')
    p_audit.add_argument('--task-id', required=True, help='Task identifier')
    p_audit.add_argument('--target-identifier', default='SPECIMEN-001', help='Target specimen identifier')
    p_audit.add_argument('--primary-metric', type=float, default=10.0, help='Primary measurement value')
    p_audit.add_argument('--secondary-metric', type=float, default=5.0, help='Secondary measurement value')
    p_audit.add_argument('--status-descriptor', default='NOMINAL', help='Status descriptor')
    p_audit.add_argument('--is-critical-flag', action='store_true', help='Mark as critical')
    p_audit.set_defaults(func=cmd_audit)

    # chat
    p_chat = subparsers.add_parser('chat', help='Query the air-gapped supervisory assistant')
    p_chat.add_argument('query', nargs='+', help='Query text')
    p_chat.set_defaults(func=cmd_chat)

    # verify-audit
    p_verify = subparsers.add_parser('verify-audit', help='Verify HMAC-SHA256 audit trail integrity')
    p_verify.set_defaults(func=cmd_verify_audit)

    args = parser.parse_args(argv)
    try:
        args.func(args)
        return 0
    except SystemExit as e:
        return e.code if isinstance(e.code, int) else 1
    except Exception:
        return 1


if __name__ == '__main__':
    sys.exit(main())
