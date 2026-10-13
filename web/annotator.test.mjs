import test from "node:test";
import assert from "node:assert/strict";
import { annotateVCF, classifyAF, toCSV } from "./annotator.mjs";

const header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO";
test("single allele INFO/AF and category", () => {
  const rows = annotateVCF(header + "\n1\t100\trs1\tA\tG\t.\tPASS\tAF=0.005");
  assert.equal(rows.length, 1);
  assert.equal(rows[0].AF, 0.005);
  assert.equal(rows[0].category, "rare");
  assert.equal(rows[0].AF_source, "INFO/AF");
});
test("multi-ALT AF values align with each ALT", () => {
  const rows = annotateVCF(header + "\n1\t100\t.\tA\tG,T\t.\tPASS\tAF=0.1,0.001");
  assert.deepEqual(rows.map(r => [r.ALT, r.AF]), [["G", 0.1], ["T", 0.001]]);
});
test("GT aggregates samples and counts only selected ALT", () => {
  const rows = annotateVCF(header + "\tFORMAT\tS1\tS2\n1\t1\t.\tA\tG,T\t.\tPASS\t.\tGT\t0/2\t1/2");
  assert.deepEqual(rows.map(r => r.AF), [0.25, 0.5]);
  assert.ok(rows.every(r => r.AF_source === "GT"));
});
test("missing AF stays unknown; zero INFO AF stays zero", () => {
  const rows = annotateVCF(header + "\n1\t1\t.\tA\tG\t.\tPASS\t.\n1\t2\t.\tA\tG\t.\tPASS\tAF=0");
  assert.deepEqual(rows.map(r => r.category), ["unknown", "ultra-rare"]);
  assert.equal(rows[0].AF, null);
});
test("AD is explicitly a read-depth fraction", () => {
  const row = annotateVCF(header + "\tFORMAT\tS1\n1\t1\t.\tA\tG\t.\tPASS\t.\tGT:AD\t./.:2,8")[0];
  assert.equal(row.AF, 0.8);
  assert.equal(row.AF_source, "AD read fraction");
});
test("invalid data fails without creating output", () => {
  assert.throws(() => annotateVCF("1\t1\t.\tA\tG\t.\tPASS\tAF=0.5"), /#CHROM/);
  assert.throws(() => annotateVCF(header + "\n1\t1\t.\tA\tG\t.\tPASS\tAF=1.2"), /invalid AF/);
  assert.throws(() => annotateVCF(header + "\n1\t0\t.\tA\tG\t.\tPASS\t."), /invalid POS/);
  assert.throws(() => annotateVCF(header + "\n1\t1\t.\tA\tG,T\t.\tPASS\tAF=0.5"), /cardinality/);
});
test("CSV quotes strings and prevents spreadsheet formula execution", () => {
  const csv = toCSV([{CHROM: '=1+1', POS: 1, ID: '"x"', REF: "A", ALT: "G", AF: null,
                      AF_source: "unavailable", category: "unknown"}]);
  assert.match(csv, /"'=1\+1"/);
  assert.match(csv, /"""x"""/);
});
test("empty VCF outputs a header-only CSV", () => {
  assert.deepEqual(annotateVCF(header), []);
  assert.equal(toCSV([]).split("\r\n")[0], "CHROM,POS,ID,REF,ALT,AF,AF_source,category");
});
test("frequency category thresholds", () => {
  assert.equal(classifyAF(0.01), "common");
  assert.equal(classifyAF(0.001), "rare");
  assert.equal(classifyAF(0.0001), "ultra-rare");
});
