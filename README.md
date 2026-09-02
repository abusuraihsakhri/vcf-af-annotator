# VCF AF Annotator

> **Domain:** Clinical Decision Support & Biomedical Computing  
> **Reference Guidelines & Standards:** `Standard Clinical Formulations & ISO/IEC Quality Frameworks`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

VCF Allele Frequency Annotator — Real Implementation

Parses VCF files, extracts/calculates allele frequencies, computes population
genetics metrics (Hardy-Weinberg, chi-square), and classifies variants by AF.

Zero external dependencies — Python stdlib only.

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Analytical Functions

- **`parse_vcf_line()`**: Parse a single VCF data line into a dict.

VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT SAMPLE...]
- **`parse_info()`**: Parse VCF INFO field into a dictionary.

Format: KEY1=val1;KEY2=val2;FLAG
- **`parse_genotype()`**: Parse VCF FORMAT/SAMPLE fields into a dictionary.

Handles GT, DP, AD, GQ, PL fields.
- **`parse_gt_alleles()`**: Parse genotype string like '0/1' or '1|0' into allele indices.
- **`extract_af_from_info()`**: Extract allele frequency from VCF INFO field (AF tag).

---

## 📐 Mathematical Formulation & Logic

```text
  Parses VCF files, extracts/calculates allele frequencies, computes population
  """Calculate allele frequency from genotype (GT) fields across samples.
  """Calculate allele frequency from AD (allelic depth) field."""
  """Calculate expected Hardy-Weinberg genotype frequencies.
  If af is not provided, it is calculated from observed genotype counts.
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --input <value> --output <value> --hom-ref <value> --het <value>
```

### Parameter Reference
- `--input`: Specifies input measurement or parameter value.
- `--output`: Specifies input measurement or parameter value.
- `--hom-ref`: Specifies input measurement or parameter value.
- `--het`: Specifies input measurement or parameter value.
- `--hom-alt`: Specifies input measurement or parameter value.
- `--af`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `query` | Parameter / observation metric | Required |
| `name` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t vcf-af-annotator .
docker run -p 8000:8000 vcf-af-annotator
```
