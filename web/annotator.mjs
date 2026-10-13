/* Local VCF allele-frequency annotation. No uploads or external dependencies. */
export const COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "AF", "AF_source", "category"];

export function classifyAF(af) {
  if (af === null) return "unknown";
  if (!Number.isFinite(af) || af < 0 || af > 1) throw new Error("AF must be in [0, 1]");
  return af >= 0.01 ? "common" : af >= 0.001 ? "rare" : "ultra-rare";
}

function numberAF(value, line) {
  if (value === undefined || value === "" || value === ".") return null;
  const af = Number(value);
  if (!Number.isFinite(af) || af < 0 || af > 1) {
    throw new Error("Line " + line + ": invalid AF '" + value + "'");
  }
  return af;
}

function fromGT(fields, altIndex) {
  if (fields.length < 10) return null;
  const format = fields[8].split(":");
  const index = format.indexOf("GT");
  if (index < 0) return null;
  let total = 0, selected = 0;
  for (const sample of fields.slice(9)) {
    const gt = (sample.split(":")[index] || "").split(/[\/|]/);
    for (const allele of gt) {
      if (allele === "." || allele === "") continue;
      if (!/^\d+$/.test(allele)) continue;
      total++;
      if (Number(allele) === altIndex) selected++;
    }
  }
  return total ? selected / total : null;
}

function fromAD(fields, altIndex) {
  if (fields.length < 10) return null;
  const index = fields[8].split(":").indexOf("AD");
  if (index < 0) return null;
  const raw = fields[9].split(":")[index];
  if (!raw) return null;
  const parts = raw.split(",");
  if (altIndex >= parts.length || parts.some(v => !/^\d+$/.test(v))) return null;
  const depths = parts.map(Number), total = depths.reduce((a, b) => a + b, 0);
  return total > 0 ? depths[altIndex] / total : null;
}

export function annotateVCF(input) {
  if (typeof input !== "string") throw new TypeError("VCF text required");
  const records = [];
  let headerSeen = false, lineNo = 0;
  for (const raw of input.replace(/^\uFEFF/, "").split(/\r?\n/)) {
    lineNo++;
    const line = raw.trim();
    if (!line) continue;
    if (line.startsWith("##")) continue;
    if (line.startsWith("#CHROM\t")) {
      headerSeen = true;
      continue;
    }
    if (line.startsWith("#")) continue;
    if (!headerSeen) throw new Error("Missing VCF #CHROM header before line " + lineNo);
    const fields = line.split("\t");
    if (fields.length < 8) throw new Error("Line " + lineNo + ": expected at least 8 tab-separated columns");
    const pos = Number(fields[1]);
    if (!Number.isSafeInteger(pos) || pos < 1) {
      throw new Error("Line " + lineNo + ": invalid POS");
    }
    const alts = fields[4].split(",");
    if (!fields[3] || alts.some(x => !x || x === ".")) {
      throw new Error("Line " + lineNo + ": invalid REF/ALT");
    }
    const infoAF = fields[7].split(";").find(part => part.startsWith("AF="));
    const afValues = infoAF === undefined ? null : infoAF.slice(3).split(",");
    if (afValues && afValues.length !== alts.length) {
      throw new Error("Line " + lineNo + ": INFO/AF cardinality must match ALT count");
    }
    for (let i = 0; i < alts.length; i++) {
      let af = afValues ? numberAF(afValues[i], lineNo) : null, source = "unavailable";
      if (af !== null) source = "INFO/AF";
      if (af === null) {
        af = fromGT(fields, i + 1);
        if (af !== null) source = "GT";
      }
      if (af === null) {
        af = fromAD(fields, i + 1);
        if (af !== null) source = "AD read fraction";
      }
      records.push({
        CHROM: fields[0], POS: pos, ID: fields[2],
        REF: fields[3], ALT: alts[i],
        AF: af, AF_source: source, category: classifyAF(af),
      });
    }
  }
  if (!headerSeen) throw new Error("Missing VCF #CHROM header");
  return records;
}

function safeCSVCell(value) {
  if (value === null || value === undefined) return "";
  let str = String(value);
  // Spreadsheet formula injection protection for user-supplied identifiers.
  if (/^[=+@\-\t\r]/.test(str)) str = "'" + str;
  return '"' + str.replace(/"/g, '""') + '"';
}

export function toCSV(rows) {
  return COLUMNS.join(",") + "\r\n" +
    rows.map(row => COLUMNS.map(key => safeCSVCell(row[key])).join(",")).join("\r\n") +
    (rows.length ? "\r\n" : "");
}
