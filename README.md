# VCF AF Annotator

VCF allele-frequency extraction and annotation for exploratory genomic analysis. A Python CLI processes VCF files and CSV tables; a static browser application annotates plain-text VCF files locally and exports CSV.

## Browser application

Open [web/index.html](web/index.html) locally, or run a static web server from the repository root and visit /web/. The repository also includes a root \`index.html\` that forwards to the browser app for compatible GitHub Pages deployments. **GitHub Pages live availability has not yet been independently verified.**

Upload or paste a plain-text VCF and select **Annotate variants**. Results can be reviewed in the table (first 100 rows) or downloaded in full. No data is uploaded by the browser app. File handling is limited to 10 MB and uncompressed VCF; the Python CLI also handles \`.vcf.gz\`.

## Python CLI

Requires Python 3.10+. Core VCF annotation uses only the standard library.

\`\`\`bash
git clone https://github.com/abusuraihsakhri/vcf-af-annotator.git
cd vcf-af-annotator
python cli.py annotate -i variants.vcf -o annotations.csv
python cli.py annotate -i variants.vcf.gz -o annotations.csv
python cli.py batch -i sample.csv -o batch.csv
python cli.py classify 0.005
python cli.py hwe --hom-ref 25 --het 50 --hom-alt 25
\`\`\`

The output CSV includes CHROM, POS, REF, a single ALT allele, AF, \`af_source\`, category and, where appropriate, HWE results. Multi-ALT VCF records yield one CSV row per ALT.

**AF priority:** (1) INFO/AF for the corresponding ALT; (2) frequency among called sample genotype alleles (GT); (3) fraction of reads supporting ALT in the first sample's AD field. AD-based fractions are labelled as such and **are not estimates of population allele frequency**. When none is present, AF remains empty and category is \`unknown\` (not ultra-rare). INFO/AF and GT may encode different biological quantities; choose input data and interpretations accordingly.

The repository's thresholds are:
- **Common:** AF ≥ 0.01
- **Rare:** 0.001 ≤ AF < 0.01
- **Ultra-rare:** AF < 0.001

These are descriptive frequency bins, **not** clinical pathogenicity classifications. HWE chi-square output is omitted for single-sample, multiallelic, incomplete/non-diploid or underpowered genotype groups. The chi-square approximation requires adequate expected counts; use an exact HWE method for small samples.

## Tests

\`\`\`bash
python -m pip install -r requirements.txt httpx
python -m pytest -q
python -m compileall -q .
node --test web/annotator.test.mjs
\`\`\`

Node.js is needed only for JavaScript regression tests; the browser application itself has no build step or external runtime dependencies.

## Optional API and container

Run a FastAPI service for the separate, rule-based supervisory/audit subsystem:

\`\`\`bash
python -m pip install -r requirements.txt
export AUDIT_SECRET_KEY="$(python -c 'import secrets; print(secrets.token_urlsafe(48))')"
export API_TOKEN="$(python -c 'import secrets; print(secrets.token_urlsafe(48))')"
uvicorn agents.api:app --host 127.0.0.1 --port 8000
\`\`\`

Endpoints: \`GET /\` (browser app redirect), \`GET /web/\` (browser app), \`GET /health\`, \`GET /metrics\`, \`POST /api/audit\`, \`POST /api/chat\`, and \`GET /api/audit/logs\`. The three \`/api/*\` routes require an \`X-API-Key\` header matching the \`API_TOKEN\` environment variable. The audit HMAC chain uses \`AUDIT_SECRET_KEY\`. Without those variables, secure endpoints are unavailable. Do not expose this demonstration API directly to untrusted networks; use TLS, access controls and a production storage model.

Docker and Docker Compose are optional; they run the same FastAPI server. Set \`API_TOKEN\` and \`AUDIT_SECRET_KEY\` in your own shell or an ignored local \`.env\` file before running \`docker compose up --build\`. Compose binds the service to localhost by default.

The supervisory subsystem uses fixed, generic measurement thresholds unrelated to VCF pathogenicity, and its chat adapter returns deterministic mock text. Its current audit trail is an **in-memory, keyed, tamper-evident chain**, not a durable external audit store. Regex identifier screening is a partial safeguard, **not guaranteed PHI de-identification or regulatory compliance**. Do not process identifiable clinical data without independent security and privacy review.

## Data handling and limitations

Browser files remain on the user's device; no analytics, file uploads or remote API calls are implemented in the browser interface. The Python CLI reads and writes local files. If the separate FastAPI server is used, supplied payloads are processed in server memory, with hashed metadata included in its audit chain. There is no automatic gnomAD download, variant consequence annotation, or clinical interpretation.

## Technology and compatibility

Python 3.10–3.12 (standard-library core), FastAPI/Pydantic for the optional API, vanilla JavaScript ES modules for recent Chrome, Edge, Firefox and Safari. The browser interface is static and can be served on GitHub Pages; Python itself is **not executed in the browser**. HWE/statistical calculations are available via the Python CLI.

## License

[MIT](LICENSE).
