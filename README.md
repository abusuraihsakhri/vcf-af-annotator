# VCF AF Annotator

> **Domain:** Clinical Decision Support & Biomedical Computing
> **Reference Guidelines & Standards:** CAP / CLSI / ISO Standards

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%20%3.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## What It Does

VCF Allele Frequency Annotator parses VCF (Variant Call Format) files, extracts and calculates allele frequencies, computes population genetics metrics (Hardy-Weinberg Equilibrium, chi-square tests), and classifies variants by allele frequency.

The core module (`vcf_af.py`) uses zero external dependencies — Python stdlib only. The enterprise features (FastAPI server, audit trail, supervisor) require additional dependencies listed in `requirements.txt`.

---

## Installation

```bash
# Clone the repository
git clone https://github.com/abusuraihsakhri/vcf-af-annotator.git
cd vcf-af-annotator

# Optional: Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies (for enterprise features)
pip install -r requirements.txt
```

---

## CLI Usage

### Annotate a VCF File
```bash
python cli.py annotate -i input.vcf -o annotated.csv
```

### Run Hardy-Weinberg Equilibrium Test
```bash
python cli.py hwe --hom-ref 25 --het 50 --hom-alt 25 --af 0.5
```

### Classify a Variant by Allele Frequency
```bash
python cli.py classify 0.05
```

### Parse a Single VCF Line
```bash
python cli.py parse "chr1\t100\t.\tA\tG\t50\tPASS\tAF=0.25\tGT:DP\t0/1:30"
```

### Batch Process from CSV
```bash
python cli.py batch -i variants.csv -o results.csv
```

### Enterprise Audit (Supervisor Dispatch)
```bash
# Set audit key (required)
export AUDIT_SECRET_KEY="your-secret-key-here"

python cli.py audit --task-id "TASK-001" --primary-metric 12.0
```

### Query the Supervisory Assistant
```bash
python cli.py chat "What is the clinical significance?"
```

### Verify Audit Trail Integrity
```bash
python cli.py verify-audit
```

---

## Key Capabilities & Algorithmic Modules

### Analytical Functions

- **`parse_vcf_line()`**: Parse a single VCF data line into a dict.
  VCF columns: CHROM POS ID REF ALT QUAL FILTER INFO [FORMAT SAMPLE...]
- **`parse_info()`**: Parse VCF INFO field into a dictionary.
  Format: KEY1=val1;KEY2=val2;FLAG
- **`parse_genotype()`**: Parse VCF FORMAT/SAMPLE fields into a dictionary.
  Handles GT, DP, AD, GQ, PL fields.
- **`parse_gt_alleles()`**: Parse genotype string like '0/1' or '1|0' into allele indices.
- **`extract_af_from_info()`**: Extract allele frequency from VCF INFO field (AF tag).
- **`calculate_af_from_genotypes()`**: Calculate allele frequency from genotype (GT) fields across samples.
- **`calculate_af_from_ad()`**: Calculate allele frequency from AD (allelic depth) field.
- **`hardy_weinberg_expected()`**: Calculate expected Hardy-Weinberg genotype frequencies.
- **`chi_square_hwe()`**: Chi-square test for Hardy-Weinberg Equilibrium deviation.
- **`classify_variant_by_af()`**: Classify variants as common, rare, or ultra-rare.

---

## Mathematical Formulation

Given allele frequency p (ref) and q (alt):
- p² — expected homozygous reference frequency
- 2pq — expected heterozygous frequency
- q² — expected homozygous alternate frequency

Variant Classification:
- Common: AF >= 1% (polymorphism)
- Rare: 0.1% <= AF < 1%
- Ultra-rare: AF < 0.1%

---

## Security Architecture

- **Zero-PHI Outbound Interceptor:** Active regex inspection blocking SSNs, MRNs, phone numbers, emails, DOBs, and patient names.
- **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation. Requires `AUDIT_SECRET_KEY` environment variable — no hardcoded defaults.
- **Path Traversal Protection:** All file I/O operations validate paths against directory traversal and null byte injection.
- **Air-Gapped LLM Reasoning Adapter:** Deterministic mock by default; pluggable for Ollama, Claude, OpenAI.
- **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## Testing

```bash
# Run all tests
pytest -v

# Run with coverage
pytest -v --cov=. --cov-report=html

# Run specific test module
pytest test_vcf_af.py -v
pytest tests/ -v
```

---

## API Server

```bash
# Start the FastAPI server
uvicorn agents.api:app --host 0.0.0.0 --port 8000

# Endpoints:
# GET  /health        - Health check
# GET  /metrics       - Prometheus metrics
# POST /api/audit     - Dispatch audit task
# POST /api/chat      - Query supervisor
# GET  /api/audit/logs - Retrieve audit trail
```

---

## Container Deployment

```bash
docker build -t vcf-af-annotator .
docker run -p 8000:8000 -e AUDIT_SECRET_KEY=your-secret-key vcf-af-annotator
```

Or with docker-compose:

```bash
docker-compose up -d
```

---

## Project Structure

```
vcf-af-annotator/
├── vcf_af.py           # Core VCF parsing & AF calculation (stdlib only)
├── cli.py              # CLI interface
├── test_vcf_af.py      # Core module tests
├── conftest.py         # Pytest configuration
├── requirements.txt    # Dependencies
├── agents/             # Enterprise supervisor system
│   ├── api.py          # FastAPI server
│   ├── base.py         # PHI guard, audit trail, security
│   ├── models.py       # Pydantic schemas
│   ├── supervisor.py   # Master coordinator
│   ├── workers.py      # Specialized domain workers
│   ├── llm_factory.py  # LLM integration
│   ├── metrics.py      # Prometheus metrics
│   ├── learning.py     # Bayesian calibration
│   └── streamer.py     # WebSocket telemetry
├── tests/              # Integration tests
│   ├── test_enrichment.py
│   └── test_vcf_af_annotator.py
├── web/                # Operations console (HTML)
├── enrichment.py       # Feature enrichment engine
├── simulator.py        # High-throughput stress tester
└── Dockerfile          # Container build
```

---

## License

MIT License. See [LICENSE](LICENSE) for details.
