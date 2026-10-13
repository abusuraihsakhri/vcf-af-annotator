"""Regression tests for AF provenance, parsing, audit integrity and API authorization."""
import csv

import pytest
from fastapi.testclient import TestClient
from vcf_af import (
    annotate_vcf_record, annotate_vcf_file, calculate_af_from_genotypes,
    extract_af_from_info, parse_info, parse_vcf_line, process_batch,
)
from agents.api import app
from agents.base import AuditTrail

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def test_multiallelic_info_per_alt(tmp_path):
    path = tmp_path / "multi.vcf"
    out = tmp_path / "result.csv"
    path.write_text(VCF_HEADER + "1\t15\t.\tA\tC,G\t.\tPASS\tAF=0.2,0.001\n")
    rows = annotate_vcf_file(str(path), str(out))
    assert [(x["ALT"], x["AF"]) for x in rows] == [("C", 0.2), ("G", 0.001)]
    assert all(x["af_source"] == "INFO/AF" for x in rows)
    with out.open(newline="") as f:
        assert len(list(csv.DictReader(f))) == 2


def test_multi_sample_genotypes_ignore_missing():
    record = parse_vcf_line("2\t11\t.\tA\tG,T\t.\tPASS\t.\tGT\t0/2\t1/2\t./.")
    first = annotate_vcf_record(record, 1)
    second = annotate_vcf_record(record, 2)
    assert first["AF"] == 0.25 and second["AF"] == 0.5
    assert first["af_source"] == "GT"
    assert "hwe_p_value" not in first


def test_ad_fallback_labeled_depth_fraction():
    rec = parse_vcf_line("1\t1\t.\tA\tC\t.\tPASS\t.\tGT:AD\t./.:1,3")
    annotation = annotate_vcf_record(rec)
    assert annotation["AF"] == 0.75
    assert annotation["af_source"] == "AD read fraction"


def test_missing_frequency_not_zero():
    rec = parse_vcf_line("1\t1\t.\tA\tC\t.\tPASS\t.")
    result = annotate_vcf_record(rec)
    assert result["AF"] is None and result["category"] == "unknown"


@pytest.mark.parametrize("af", ["-0.1", "1.1", "nan", "inf"])
def test_invalid_af_rejected(af):
    rec = parse_vcf_line(f"1\t1\t.\tA\tG\t.\tPASS\tAF={af}")
    with pytest.raises(ValueError):
        annotate_vcf_record(rec)


def test_empty_vcf_writes_header(tmp_path):
    input_file, output_file = tmp_path / "empty.vcf", tmp_path / "out.csv"
    input_file.write_text(VCF_HEADER)
    assert annotate_vcf_file(str(input_file), str(output_file)) == []
    assert output_file.read_text().startswith("CHROM,POS,REF,ALT,AF,")


def test_gzip_vcf(tmp_path):
    import gzip
    input_file, output_file = tmp_path / "in.vcf.gz", tmp_path / "out.csv"
    with gzip.open(input_file, "wt") as f:
        f.write(VCF_HEADER + "1\t1\t.\tA\tC\t.\tPASS\tAF=0.15\n")
    assert annotate_vcf_file(str(input_file), str(output_file))[0]["AF"] == 0.15


def test_missing_gt_not_double_counted():
    records = [{"GENOTYPE": {"alleles": [-1, -1]}},
               {"GENOTYPE": {"alleles": [0, 1]}}]
    res = calculate_af_from_genotypes(records)
    assert res["total_alleles"] == 2
    assert res["genotype_counts"]["missing"] == 1


def test_tamper_evident_audit_recomputes_signature():
    audit = AuditTrail("private-unit-test-key")
    audit.log("worker", "system", "event", {"status": "OK"})
    assert audit.verify_integrity()
    exported = audit.get_trail()
    exported[0]["actor"] = "modified"
    assert audit.verify_integrity()  # Returned entries cannot mutate internal logs
    audit.logs[0]["actor"] = "modified"
    assert not audit.verify_integrity()


def test_api_requires_token(monkeypatch):
    client = TestClient(app)
    monkeypatch.delenv("API_TOKEN", raising=False)
    payload = {"task_id": "T1", "target_identifier": "KEY-01", "primary_metric": 1.0}
    assert client.post("/api/audit", json=payload).status_code == 503
    monkeypatch.setenv("API_TOKEN", "test-only-secret")
    assert client.post("/api/audit", json=payload).status_code == 401
    assert client.get("/api/audit/logs").status_code == 401
    result = client.post("/api/audit", json=payload,
                         headers={"X-API-Key": "test-only-secret"})
    assert result.status_code == 200
    assert result.json()["task_id"] == "T1"


def test_api_rejects_nonfinite_metric(monkeypatch):
    client = TestClient(app)
    monkeypatch.setenv("API_TOKEN", "test-only-secret")
    response = client.post("/api/audit",
        json={"task_id": "T", "target_identifier": "KEY", "primary_metric": "NaN"},
        headers={"X-API-Key": "test-only-secret"})
    assert response.status_code == 422
