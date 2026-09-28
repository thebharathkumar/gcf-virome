"""Validate the REDCap and SRA metadata deliverables.

Three classes of problem are checked, which are the three that actually cause
rejected submissions and unusable study databases:

  1. missing required fields
  2. sample IDs that disagree across files
  3. values that are invalid for their declared type or vocabulary

Exit status is 0 when there are no errors and 1 otherwise, so the workflow
fails loudly rather than producing a report nobody reads. Warnings do not
fail the run.

Field names and formats are those recorded in docs/metadata_standards.md,
which cites the REDCap and NCBI documentation they came from.
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

# --------------------------------------------------------------------------
# Schema constants. See docs/metadata_standards.md for the source of each.
# --------------------------------------------------------------------------

REDCAP_DICT_COLUMNS = [
    "Variable / Field Name",
    "Form Name",
    "Section Header",
    "Field Type",
    "Field Label",
    "Choices, Calculations, OR Slider Labels",
    "Field Note",
    "Text Validation Type OR Show Slider Number",
    "Text Validation Min",
    "Text Validation Max",
    "Identifier?",
    "Branching Logic (Show field only if...)",
    "Required Field?",
    "Custom Alignment",
    "Question Number (surveys only)",
    "Matrix Group Name",
    "Matrix Ranking?",
    "Field Annotation",
]

REDCAP_FIELD_TYPES = {
    "text", "notes", "dropdown", "radio", "checkbox", "file", "calc",
    "sql", "descriptive", "slider", "yesno", "truefalse",
}

REDCAP_VALIDATION_TYPES = {
    "", "date_ymd", "date_mdy", "date_dmy",
    "datetime_ymd", "datetime_mdy", "datetime_dmy",
    "datetime_seconds_ymd", "datetime_seconds_mdy", "datetime_seconds_dmy",
    "time", "email", "integer", "number", "phone", "zipcode", "signature",
    "number_1dp", "number_2dp", "number_3dp", "number_4dp",
}

# REDCap variable names: lowercase letters, digits and underscores, starting
# with a letter, minimum 2 characters, 26 characters recommended maximum.
REDCAP_VARNAME = re.compile(r"^[a-z][a-z0-9_]{1,25}$")

# INSDC / NCBI accepted missing-value vocabulary.
INSDC_MISSING = {
    "not applicable", "not collected", "not provided", "missing",
    "restricted access",
}

# BioSample package MIMS.me.human-oral.6.0. The seven package-mandatory
# attributes plus the two base columns that every package carries.
BIOSAMPLE_REQUIRED_COLUMNS = [
    "sample_name",
    "organism",
    "collection_date",
    "env_broad_scale",
    "env_local_scale",
    "env_medium",
    "geo_loc_name",
    "host",
    "lat_lon",
]

# SRA library/run sheet, column names taken from NCBI's own SRA metadata
# spreadsheet.
SRA_REQUIRED_COLUMNS = [
    "bioproject_accession",
    "biosample_accession",
    "library_ID",
    "title",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "platform",
    "instrument_model",
    "filetype",
]

SRA_PLATFORMS = {
    "_LS454", "ILLUMINA", "HELICOS", "ABI_SOLID", "COMPLETE_GENOMICS",
    "PACBIO_SMRT", "ION_TORRENT", "CAPILLARY", "OXFORD_NANOPORE",
}
SRA_FILETYPES = {
    "bam", "srf", "sff", "fastq", "454_native", "Helicos_native",
    "SOLiD_native", "PacBio_HDF5", "CompleteGenomics_native",
}

# An ENVO/UBERON style term carries its ID in brackets.
ONTOLOGY_TERM = re.compile(r"^.+\[[A-Za-z]+:\d+\]$")

SRA_LIBRARY_STRATEGY = {
    "WGS", "WGA", "WXS", "RNA-Seq", "miRNA-Seq", "AMPLICON", "CLONE",
    "POOLCLONE", "CLONEEND", "FINISHING", "ChIP-Seq", "MNase-Seq", "DNase-Hypersensitivity",
    "Bisulfite-Seq", "EST", "FL-cDNA", "CTS", "MRE-Seq", "MeDIP-Seq", "MBD-Seq",
    "Tn-Seq", "VALIDATION", "FAIRE-seq", "SELEX", "RIP-Seq", "ChIA-PET",
    "Synthetic-Long-Read", "Targeted-Capture", "Tethered Chromatin Conformation Capture",
    "OTHER",
}
SRA_LIBRARY_SOURCE = {
    "GENOMIC", "TRANSCRIPTOMIC", "METAGENOMIC", "METATRANSCRIPTOMIC",
    "SYNTHETIC", "VIRAL RNA", "GENOMIC SINGLE CELL",
    "TRANSCRIPTOMIC SINGLE CELL", "OTHER",
}
SRA_LIBRARY_SELECTION = {
    "RANDOM", "PCR", "RANDOM PCR", "RT-PCR", "HMPR", "MF", "CF-S", "CF-M",
    "CF-H", "CF-T", "MDA", "MSLL", "cDNA", "ChIP", "MNase", "DNAse",
    "Hybrid Selection", "Reduced Representation", "Restriction Digest",
    "5-methylcytidine antibody", "MBD2 protein methyl-CpG binding domain",
    "CAGE", "RACE", "size fractionation", "Padlock probes capture method",
    "other", "unspecified", "cDNA_oligo_dT", "cDNA_randomPriming",
    "Inverse rRNA", "Oligo-dT", "PolyA", "repeat fractionation",
}
SRA_LIBRARY_LAYOUT = {"single", "paired"}

DATE_YMD = re.compile(r"^\d{4}-\d{2}-\d{2}$")
# NCBI accepts a year, a year-month, a full date, a range, or ISO 8601.
NCBI_DATE = re.compile(
    r"^(\d{4}(-\d{2}(-\d{2})?)?|\d{4}(-\d{2}(-\d{2})?)?/\d{4}(-\d{2}(-\d{2})?)?)$"
)
# "Country" or "Country:Region"
GEO_LOC = re.compile(r"^[A-Za-z][A-Za-z .'\-]*(:\s?.+)?$")
# Decimal degrees with hemisphere letters, e.g. "50.8802 N 4.6934 E"
LAT_LON = re.compile(r"^\d+(\.\d+)?\s+[NS]\s+\d+(\.\d+)?\s+[EW]$")


class Report:
    def __init__(self):
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.notes: list[str] = []
        self.checks_run = 0

    def error(self, where, msg):
        self.errors.append(f"[{where}] {msg}")

    def warn(self, where, msg):
        self.warnings.append(f"[{where}] {msg}")

    def note(self, msg):
        self.notes.append(msg)

    def check(self):
        self.checks_run += 1


def load_csv(path, delimiter=","):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        rows = list(reader)
        return reader.fieldnames or [], rows


def is_blank(value):
    return value is None or str(value).strip() == ""


# --------------------------------------------------------------------------
# 1. REDCap data dictionary
# --------------------------------------------------------------------------

def parse_choices(raw):
    """'1, Male | 2, Female' -> {'1': 'Male', '2': 'Female'}."""
    out = {}
    for chunk in str(raw).split("|"):
        chunk = chunk.strip()
        if not chunk:
            continue
        code, sep, label = chunk.partition(",")
        if not sep:
            return None
        out[code.strip()] = label.strip()
    return out or None


def validate_dictionary(path, rep):
    cols, rows = load_csv(path)
    rep.check()
    if cols != REDCAP_DICT_COLUMNS:
        missing = [c for c in REDCAP_DICT_COLUMNS if c not in cols]
        extra = [c for c in cols if c not in REDCAP_DICT_COLUMNS]
        if missing or extra:
            rep.error("dictionary", f"column mismatch. missing={missing} unexpected={extra}")
        else:
            rep.error("dictionary", "columns present but in the wrong order")

    if not rows:
        rep.error("dictionary", "no fields defined")
        return {}

    fields = {}
    seen = set()
    for i, row in enumerate(rows, start=2):
        name = (row.get("Variable / Field Name") or "").strip()
        rep.check()
        if not name:
            rep.error("dictionary", f"line {i}: empty variable name")
            continue
        if not REDCAP_VARNAME.match(name):
            rep.error(
                "dictionary",
                f"line {i}: '{name}' is not a valid REDCap variable name "
                "(lowercase, start with a letter, <= 26 chars)",
            )
        if name in seen:
            rep.error("dictionary", f"line {i}: duplicate variable name '{name}'")
        seen.add(name)

        ftype = (row.get("Field Type") or "").strip()
        if ftype not in REDCAP_FIELD_TYPES:
            rep.error("dictionary", f"'{name}': invalid Field Type '{ftype}'")

        if is_blank(row.get("Field Label")):
            rep.error("dictionary", f"'{name}': Field Label is required")

        if is_blank(row.get("Form Name")):
            rep.error("dictionary", f"'{name}': Form Name is required")

        vtype = (row.get("Text Validation Type OR Show Slider Number") or "").strip()
        if vtype not in REDCAP_VALIDATION_TYPES:
            rep.error("dictionary", f"'{name}': unknown validation type '{vtype}'")

        choices = None
        raw_choices = row.get("Choices, Calculations, OR Slider Labels") or ""
        if ftype in {"radio", "dropdown", "checkbox"}:
            choices = parse_choices(raw_choices)
            if choices is None:
                rep.error(
                    "dictionary",
                    f"'{name}': {ftype} field needs choices as "
                    "'code, label | code, label'",
                )
        elif not is_blank(raw_choices) and ftype not in {"calc", "slider", "sql"}:
            rep.warn("dictionary", f"'{name}': choices given for field type '{ftype}'")

        required = (row.get("Required Field?") or "").strip().lower()
        if required not in {"", "y"}:
            rep.error(
                "dictionary",
                f"'{name}': Required Field? must be 'y' or empty, found '{required}'",
            )

        identifier = (row.get("Identifier?") or "").strip().lower()
        if identifier not in {"", "y"}:
            rep.error(
                "dictionary",
                f"'{name}': Identifier? must be 'y' or empty, found '{identifier}'",
            )

        fields[name] = {
            "type": ftype,
            "validation": vtype,
            "choices": choices,
            "required": required == "y",
            "min": (row.get("Text Validation Min") or "").strip(),
            "max": (row.get("Text Validation Max") or "").strip(),
            "order": i,
        }

    rep.note(f"dictionary: {len(fields)} fields across "
             f"{len({r.get('Form Name') for r in rows})} form(s)")
    return fields


# --------------------------------------------------------------------------
# 2. REDCap records
# --------------------------------------------------------------------------

def validate_records(path, fields, rep):
    cols, rows = load_csv(path)
    rep.check()
    if not cols:
        rep.error("records", "file has no header")
        return []

    record_id_field = next(iter(fields)) if fields else None
    if record_id_field and cols[0] != record_id_field:
        rep.error(
            "records",
            f"first column must be the record ID field '{record_id_field}', "
            f"found '{cols[0]}'",
        )

    for col in cols:
        rep.check()
        base = col.split("___")[0]
        if base not in fields:
            rep.error("records", f"column '{col}' is not defined in the data dictionary")
        elif "___" in col and fields[base]["type"] != "checkbox":
            rep.error(
                "records",
                f"column '{col}' uses checkbox ___ syntax but '{base}' is "
                f"type '{fields[base]['type']}'",
            )

    if not rows:
        rep.error("records", "no data rows")
        return []

    ids = []
    for i, row in enumerate(rows, start=2):
        rid = (row.get(record_id_field) or "").strip() if record_id_field else ""
        ids.append(rid)
        rep.check()
        if not rid:
            rep.error("records", f"line {i}: empty record ID")

        for col, value in row.items():
            if col is None or col not in cols:
                continue
            base = col.split("___")[0]
            spec = fields.get(base)
            if spec is None:
                continue

            if is_blank(value):
                if spec["required"]:
                    rep.error("records", f"line {i} ({rid}): required field '{col}' is empty")
                continue

            value = str(value).strip()

            if spec["type"] in {"radio", "dropdown"} and spec["choices"]:
                if value not in spec["choices"]:
                    rep.error(
                        "records",
                        f"line {i} ({rid}): '{col}' = '{value}' is not a valid code "
                        f"(expected one of {sorted(spec['choices'])})",
                    )
            elif spec["type"] == "checkbox":
                if value not in {"0", "1"}:
                    rep.error(
                        "records",
                        f"line {i} ({rid}): checkbox '{col}' must be 0 or 1, found '{value}'",
                    )
            elif spec["type"] == "yesno":
                if value not in {"0", "1"}:
                    rep.error(
                        "records",
                        f"line {i} ({rid}): yesno '{col}' must be 0 or 1, found '{value}'",
                    )

            vt = spec["validation"]
            if vt == "date_ymd" and not DATE_YMD.match(value):
                rep.error(
                    "records",
                    f"line {i} ({rid}): '{col}' = '{value}' is not YYYY-MM-DD",
                )
            elif vt == "integer":
                try:
                    number = int(value)
                except ValueError:
                    rep.error("records", f"line {i} ({rid}): '{col}' = '{value}' is not an integer")
                else:
                    _range_check(rep, i, rid, col, number, spec)
            elif vt.startswith("number"):
                try:
                    number = float(value)
                except ValueError:
                    rep.error("records", f"line {i} ({rid}): '{col}' = '{value}' is not a number")
                else:
                    _range_check(rep, i, rid, col, number, spec)

    dupes = {v for v in ids if ids.count(v) > 1 and v}
    if dupes:
        rep.error("records", f"duplicate record IDs: {sorted(dupes)}")

    rep.note(f"records: {len(rows)} rows, {len(cols)} columns")
    return ids


def _range_check(rep, line, rid, col, number, spec):
    for bound, comp, word in ((spec["min"], float.__lt__, "below"),
                              (spec["max"], float.__gt__, "above")):
        if bound:
            try:
                limit = float(bound)
            except ValueError:
                continue
            if comp(float(number), limit):
                rep.error(
                    "records",
                    f"line {line} ({rid}): '{col}' = {number} is {word} the declared "
                    f"limit {bound}",
                )


# --------------------------------------------------------------------------
# 3. SRA / BioSample sheet
# --------------------------------------------------------------------------

def validate_biosample(path, rep):
    """BioSample attributes sheet, package MIMS.me.human-oral.6.0."""
    cols, rows = load_csv(path, delimiter="\t")
    rep.check()

    missing = [c for c in BIOSAMPLE_REQUIRED_COLUMNS if c not in cols]
    if missing:
        rep.error("biosample", f"missing required columns: {missing}")

    if not rows:
        rep.error("biosample", "no data rows")
        return []

    names = []
    for i, row in enumerate(rows, start=2):
        name = (row.get("sample_name") or "").strip()
        names.append(name)
        rep.check()
        if not name:
            rep.error("biosample", f"line {i}: empty sample_name")

        for col in BIOSAMPLE_REQUIRED_COLUMNS:
            if col not in cols:
                continue
            if not (row.get(col) or "").strip():
                rep.error(
                    "biosample",
                    f"line {i} ({name}): mandatory attribute '{col}' is empty. NCBI "
                    f"requires a value or one of {sorted(INSDC_MISSING)}",
                )

        date = (row.get("collection_date") or "").strip()
        if date and date.lower() not in INSDC_MISSING and not NCBI_DATE.match(date):
            rep.error(
                "biosample",
                f"line {i} ({name}): collection_date '{date}' is not an accepted "
                "NCBI date or start/end range",
            )

        geo = (row.get("geo_loc_name") or "").strip()
        if geo and geo.lower() not in INSDC_MISSING and not GEO_LOC.match(geo):
            rep.error(
                "biosample",
                f"line {i} ({name}): geo_loc_name '{geo}' is not 'Country' "
                "or 'Country:Region'",
            )

        ll = (row.get("lat_lon") or "").strip()
        if ll and ll.lower() not in INSDC_MISSING and not LAT_LON.match(ll):
            rep.error(
                "biosample",
                f"line {i} ({name}): lat_lon '{ll}' is not decimal degrees with a "
                "compass direction, for example '50.8802 N 4.6934 E'",
            )

        for col in ("env_broad_scale", "env_local_scale", "env_medium"):
            value = (row.get(col) or "").strip()
            rep.check()
            if value and value.lower() not in INSDC_MISSING and not ONTOLOGY_TERM.match(value):
                rep.warn(
                    "biosample",
                    f"line {i} ({name}): {col} '{value}' has no ontology ID in brackets; "
                    "MIxS expects a term such as 'oral cavity [UBERON:0000167]'",
                )

    dupes = {v for v in names if names.count(v) > 1 and v}
    if dupes:
        rep.error("biosample", f"duplicate sample_name values: {sorted(dupes)}")

    rep.note(f"biosample: {len(rows)} rows, {len(cols)} columns")
    return names


def validate_sra(path, rep):
    """SRA library and run sheet."""
    cols, rows = load_csv(path, delimiter="\t")
    rep.check()

    missing = [c for c in SRA_REQUIRED_COLUMNS if c not in cols]
    if missing:
        rep.error("sra", f"missing required columns: {missing}")

    if not rows:
        rep.error("sra", "no data rows")
        return []

    names = []
    for i, row in enumerate(rows, start=2):
        name = (row.get("library_ID") or "").strip()
        names.append(name)
        rep.check()
        if not name:
            rep.error("sra", f"line {i}: empty library_ID")

        for col in SRA_REQUIRED_COLUMNS:
            if col not in cols:
                continue
            if not (row.get(col) or "").strip():
                rep.error("sra", f"line {i} ({name}): required column '{col}' is empty")

        for col, allowed in (
            ("library_strategy", SRA_LIBRARY_STRATEGY),
            ("library_source", SRA_LIBRARY_SOURCE),
            ("library_selection", SRA_LIBRARY_SELECTION),
            ("library_layout", SRA_LIBRARY_LAYOUT),
            ("platform", SRA_PLATFORMS),
            ("filetype", SRA_FILETYPES),
        ):
            value = (row.get(col) or "").strip()
            rep.check()
            if value and value not in allowed:
                rep.error("sra", f"line {i} ({name}): {col} '{value}' is not an allowed term")

        for col in ("filename", "filename2"):
            value = (row.get(col) or "").strip()
            if value and not value.endswith((".fastq.gz", ".fastq", ".bam")):
                rep.warn("sra", f"line {i} ({name}): {col} '{value}' has an unusual extension")

    dupes = {v for v in names if names.count(v) > 1 and v}
    if dupes:
        rep.error("sra", f"duplicate library_ID values: {sorted(dupes)}")

    rep.note(f"sra: {len(rows)} rows, {len(cols)} columns")
    return names


# --------------------------------------------------------------------------
# 4. Cross-file agreement
# --------------------------------------------------------------------------

def validate_cross(manifest_ids, record_ids, biosample_names, sra_names, rep):
    rep.check()
    m = set(manifest_ids)
    r = set(record_ids)
    b = set(biosample_names)
    s = set(sra_names)

    for label, other in (
        ("REDCap records", r),
        ("BioSample sheet", b),
        ("SRA sheet", s),
    ):
        only_manifest = sorted(m - other)
        only_other = sorted(other - m)
        if only_manifest:
            rep.error("cross-file", f"in manifest but not in {label}: {only_manifest}")
        if only_other:
            rep.error("cross-file", f"in {label} but not in manifest: {only_other}")

    if b and s and b != s:
        rep.error(
            "cross-file",
            f"BioSample and SRA sheets disagree: "
            f"biosample-only={sorted(b - s)} sra-only={sorted(s - b)}",
        )
    if r and s and r != s:
        rep.error(
            "cross-file",
            f"REDCap records and SRA sheet disagree: "
            f"records-only={sorted(r - s)} sra-only={sorted(s - r)}",
        )

    if m and m == r == b == s:
        rep.note(
            f"cross-file: {len(m)} sample IDs agree across manifest, REDCap, "
            "BioSample and SRA"
        )


# --------------------------------------------------------------------------

def run(dictionary, records, biosample, sra, manifest, report_path=None):
    rep = Report()

    fields = validate_dictionary(dictionary, rep)
    record_ids = validate_records(records, fields, rep)
    biosample_names = validate_biosample(biosample, rep)
    sra_names = validate_sra(sra, rep)

    manifest_ids = []
    if manifest:
        _, rows = load_csv(manifest, delimiter="\t")
        manifest_ids = [r["sample_id"] for r in rows if r.get("sample_id")]

    validate_cross(manifest_ids, record_ids, biosample_names, sra_names, rep)

    lines = ["METADATA VALIDATION REPORT", "=" * 72, ""]
    lines.append(f"data dictionary : {dictionary}")
    lines.append(f"records         : {records}")
    lines.append(f"biosample sheet : {biosample}")
    lines.append(f"sra sheet       : {sra}")
    lines.append(f"manifest        : {manifest}")
    lines.append("")
    for note in rep.notes:
        lines.append(f"  {note}")
    lines.append("")
    lines.append(f"checks run : {rep.checks_run}")
    lines.append(f"errors     : {len(rep.errors)}")
    lines.append(f"warnings   : {len(rep.warnings)}")
    lines.append("")

    if rep.errors:
        lines.append("ERRORS")
        lines.append("-" * 72)
        lines.extend(f"  {e}" for e in rep.errors)
        lines.append("")
    if rep.warnings:
        lines.append("WARNINGS")
        lines.append("-" * 72)
        lines.extend(f"  {w}" for w in rep.warnings)
        lines.append("")

    lines.append("PASS" if not rep.errors else "FAIL")
    text = "\n".join(lines) + "\n"

    print(text)
    if report_path:
        Path(report_path).parent.mkdir(parents=True, exist_ok=True)
        Path(report_path).write_text(text)

    return rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dictionary", required=True)
    ap.add_argument("--records", required=True)
    ap.add_argument("--biosample", required=True)
    ap.add_argument("--sra", required=True)
    ap.add_argument("--manifest", default=None)
    ap.add_argument("--report", default=None)
    args = ap.parse_args(argv)

    rep = run(
        args.dictionary, args.records, args.biosample, args.sra,
        args.manifest, args.report,
    )
    return 1 if rep.errors else 0


if __name__ == "__main__":
    sys.exit(main())
