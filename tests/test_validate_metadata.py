"""Tests for the metadata validator.

Each test builds a small, valid metadata set, breaks exactly one thing, and
asserts that the validator catches that specific problem. A validator that
never fails is worthless, so the negative cases carry most of the weight here.
"""
import csv

import pytest

import validate_metadata as vm

# name, form, type, label, choices, validation, required, identifier, min, max
DICT_ROWS = [
    ("record_id", "participant", "text", "Sample ID", "", "", "y", "", "", ""),
    ("country", "participant", "radio", "Country", "1, Belgium | 2, Chile", "", "y", "", "", ""),
    ("diagnosis", "participant", "radio", "Diagnosis",
     "1, Periodontally healthy | 2, Periodontitis", "", "y", "", "", ""),
    ("collected_on", "specimen", "text", "Collection date", "", "date_ymd", "", "", "", ""),
    ("depth_reads", "specimen", "text", "Reads", "", "integer", "", "", "0", "1000000"),
    ("consented", "specimen", "yesno", "Consent on file", "", "", "", "", "", ""),
    ("sites", "specimen", "checkbox", "Sites sampled", "1, Molar | 2, Incisor", "", "", "", "", ""),
]

SAMPLES = ["BEHE003", "BEPE001"]

RECORD_COLUMNS = [
    "record_id", "country", "diagnosis", "collected_on",
    "depth_reads", "consented", "sites___1", "sites___2",
]
RECORD_ROWS = [
    ["BEHE003", "1", "1", "2019-03-04", "1000", "1", "1", "0"],
    ["BEPE001", "1", "2", "2019-03-05", "2000", "1", "0", "1"],
]

BIOSAMPLE_BASE = {
    "sample_name": "",
    "organism": "human oral metagenome",
    "collection_date": "2017/2019",
    "env_broad_scale": "oral cavity [UBERON:0000167]",
    "env_local_scale": "gingiva [UBERON:0001828]",
    "env_medium": "bodily fluid material [ENVO:02000019]",
    "geo_loc_name": "Belgium",
    "host": "Homo sapiens",
    "lat_lon": "50.88025 N 4.693472 E",
}

SRA_BASE = {
    "bioproject_accession": "PRJNA1183294",
    "biosample_accession": "SAMN44622345",
    "library_ID": "",
    "title": "WGS of human oral metagenome",
    "library_strategy": "WGS",
    "library_source": "METAGENOMIC",
    "library_selection": "size fractionation",
    "library_layout": "paired",
    "platform": "ILLUMINA",
    "instrument_model": "Illumina NovaSeq 6000",
    "filetype": "fastq",
    "filename": "x_1.fastq.gz",
    "filename2": "x_2.fastq.gz",
}


def write_dictionary(path, rows=None, columns=None):
    rows = DICT_ROWS if rows is None else rows
    columns = vm.REDCAP_DICT_COLUMNS if columns is None else columns
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        for name, form, ftype, label, ch, val, req, ident, vmin, vmax in rows:
            w.writerow([name, form, "", ftype, label, ch, "", val, vmin, vmax,
                        ident, "", req, "", "", "", "", ""])


def write_records(path, rows=None, columns=None):
    columns = columns or RECORD_COLUMNS
    rows = RECORD_ROWS if rows is None else rows
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(columns)
        w.writerows(rows)


def _sheet(path, columns, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def write_biosample(path, overrides=None, drop=None, ids=None, columns=None):
    columns = columns or list(vm.BIOSAMPLE_REQUIRED_COLUMNS)
    if drop:
        columns = [c for c in columns if c not in drop]
    rows = []
    for s in (ids or SAMPLES):
        r = dict(BIOSAMPLE_BASE)
        r["sample_name"] = s
        r.update(overrides or {})
        rows.append(r)
    _sheet(path, columns, rows)


def write_sra(path, overrides=None, drop=None, ids=None, columns=None):
    columns = columns or list(vm.SRA_REQUIRED_COLUMNS) + ["filename", "filename2"]
    if drop:
        columns = [c for c in columns if c not in drop]
    rows = []
    for s in (ids or SAMPLES):
        r = dict(SRA_BASE)
        r["library_ID"] = s
        r.update(overrides or {})
        rows.append(r)
    _sheet(path, columns, rows)


def write_manifest(path, ids=None):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "run_accession"])
        for i, s in enumerate(ids or SAMPLES):
            w.writerow([s, f"SRR{1000 + i}"])


@pytest.fixture
def meta(tmp_path):
    """A complete, valid metadata set. Each test breaks one piece of it."""
    paths = {
        "dictionary": tmp_path / "dict.csv",
        "records": tmp_path / "records.csv",
        "biosample": tmp_path / "biosample.tsv",
        "sra": tmp_path / "sra.tsv",
        "manifest": tmp_path / "manifest.tsv",
    }
    write_dictionary(paths["dictionary"])
    write_records(paths["records"])
    write_biosample(paths["biosample"])
    write_sra(paths["sra"])
    write_manifest(paths["manifest"])
    return paths


def run(meta):
    return vm.run(
        str(meta["dictionary"]), str(meta["records"]),
        str(meta["biosample"]), str(meta["sra"]), str(meta["manifest"]),
    )


def errors_matching(rep, needle):
    return [e for e in rep.errors if needle.lower() in e.lower()]


# ---------------------------------------------------------------- baseline

def test_valid_metadata_passes(meta):
    rep = run(meta)
    assert rep.errors == [], f"expected a clean pass, got: {rep.errors}"
    assert rep.checks_run > 0


# ------------------------------------------------- 1. missing required fields

def test_missing_required_value_in_records(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "", "1", "2019-03-04", "1000", "1", "1", "0"],
        ["BEPE001", "1", "2", "2019-03-05", "2000", "1", "0", "1"],
    ])
    rep = run(meta)
    assert errors_matching(rep, "required field 'country' is empty")


def test_missing_mandatory_biosample_column(meta):
    write_biosample(meta["biosample"], drop=["env_medium"])
    rep = run(meta)
    assert errors_matching(rep, "missing required columns")
    assert errors_matching(rep, "env_medium")


def test_empty_mandatory_biosample_value(meta):
    write_biosample(meta["biosample"], overrides={"geo_loc_name": ""})
    rep = run(meta)
    assert errors_matching(rep, "mandatory attribute 'geo_loc_name' is empty")


def test_missing_required_sra_column(meta):
    write_sra(meta["sra"], drop=["library_strategy"])
    rep = run(meta)
    assert errors_matching(rep, "missing required columns")
    assert errors_matching(rep, "library_strategy")


def test_insdc_missing_terms_are_accepted(meta):
    """'not collected' is a legitimate value, not an error."""
    write_biosample(
        meta["biosample"],
        overrides={"collection_date": "not collected", "lat_lon": "not collected"},
    )
    rep = run(meta)
    assert not errors_matching(rep, "collection_date")
    assert not errors_matching(rep, "lat_lon")


# ------------------------------------------- 2. inconsistent sample IDs

def test_sample_id_missing_from_records(meta):
    write_records(meta["records"], rows=[RECORD_ROWS[0]])
    rep = run(meta)
    assert errors_matching(rep, "in manifest but not in redcap records")
    assert errors_matching(rep, "BEPE001")


def test_sample_id_extra_in_sra(meta):
    write_manifest(meta["manifest"], ids=["BEHE003"])
    rep = run(meta)
    assert errors_matching(rep, "in sra sheet but not in manifest")


def test_biosample_and_sra_disagree(meta):
    write_biosample(meta["biosample"], ids=["BEHE003", "TYPO001"])
    rep = run(meta)
    assert errors_matching(rep, "biosample") and errors_matching(rep, "TYPO001")


def test_records_and_sra_disagree(meta):
    write_records(meta["records"], rows=[
        RECORD_ROWS[0],
        ["TYPO001", "1", "2", "2019-03-05", "2000", "1", "0", "1"],
    ])
    rep = run(meta)
    assert errors_matching(rep, "TYPO001")


def test_duplicate_record_id(meta):
    write_records(meta["records"], rows=[
        RECORD_ROWS[0],
        ["BEHE003", "1", "2", "2019-03-05", "2000", "1", "0", "1"],
    ])
    rep = run(meta)
    assert errors_matching(rep, "duplicate record ids")


def test_duplicate_biosample_name(meta):
    write_biosample(meta["biosample"], ids=["BEHE003", "BEHE003"])
    rep = run(meta)
    assert errors_matching(rep, "duplicate sample_name")


def test_duplicate_sra_library_id(meta):
    write_sra(meta["sra"], ids=["BEHE003", "BEHE003"])
    rep = run(meta)
    assert errors_matching(rep, "duplicate library_id")


# ------------------------------------------------- 3. invalid values

def test_invalid_radio_code(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "9", "1", "2019-03-04", "1000", "1", "1", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "is not a valid code")


def test_radio_label_instead_of_code_is_rejected(meta):
    """A classic import failure: the label pasted in place of the code."""
    write_records(meta["records"], rows=[
        ["BEHE003", "Belgium", "1", "2019-03-04", "1000", "1", "1", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "is not a valid code")


def test_bad_date_format(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "1", "1", "04/03/2019", "1000", "1", "1", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "is not YYYY-MM-DD")


def test_non_integer_in_integer_field(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "1", "1", "2019-03-04", "many", "1", "1", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "is not an integer")


def test_value_above_declared_maximum(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "1", "1", "2019-03-04", "99999999", "1", "1", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "above the declared limit")


def test_yesno_must_be_zero_or_one(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "1", "1", "2019-03-04", "1000", "yes", "1", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "yesno")


def test_checkbox_must_be_zero_or_one(meta):
    write_records(meta["records"], rows=[
        ["BEHE003", "1", "1", "2019-03-04", "1000", "1", "Molar", "0"],
        RECORD_ROWS[1],
    ])
    rep = run(meta)
    assert errors_matching(rep, "checkbox")


def test_undefined_column_in_records(meta):
    write_records(
        meta["records"],
        columns=RECORD_COLUMNS + ["smoking"],
        rows=[r + ["1"] for r in RECORD_ROWS],
    )
    rep = run(meta)
    assert errors_matching(rep, "not defined in the data dictionary")


def test_checkbox_syntax_on_non_checkbox_field(meta):
    cols = list(RECORD_COLUMNS)
    cols[1] = "country___1"
    write_records(meta["records"], columns=cols, rows=RECORD_ROWS)
    rep = run(meta)
    assert errors_matching(rep, "___ syntax")


def test_bad_lat_lon_format(meta):
    """Signed decimals are a common mistake; NCBI wants compass directions."""
    write_biosample(meta["biosample"], overrides={"lat_lon": "50.88, 4.69"})
    rep = run(meta)
    assert errors_matching(rep, "lat_lon")


def test_bad_collection_date_in_biosample(meta):
    write_biosample(meta["biosample"], overrides={"collection_date": "March 2019"})
    rep = run(meta)
    assert errors_matching(rep, "collection_date")


def test_invalid_library_strategy(meta):
    write_sra(meta["sra"], overrides={"library_strategy": "shotgun"})
    rep = run(meta)
    assert errors_matching(rep, "library_strategy")


def test_invalid_library_source(meta):
    write_sra(meta["sra"], overrides={"library_source": "metagenome"})
    rep = run(meta)
    assert errors_matching(rep, "library_source")


def test_invalid_platform(meta):
    write_sra(meta["sra"], overrides={"platform": "Illumina"})
    rep = run(meta)
    assert errors_matching(rep, "platform")


def test_invalid_filetype(meta):
    write_sra(meta["sra"], overrides={"filetype": "fq.gz"})
    rep = run(meta)
    assert errors_matching(rep, "filetype")


def test_ontology_term_without_id_warns(meta):
    write_biosample(meta["biosample"], overrides={"env_medium": "gingival crevicular fluid"})
    rep = run(meta)
    assert any("env_medium" in w for w in rep.warnings)
    assert not errors_matching(rep, "env_medium")


# ------------------------------------------- 4. data dictionary problems

def test_dictionary_column_mismatch(meta):
    bad = list(vm.REDCAP_DICT_COLUMNS)
    bad[10] = "Identifier"  # question mark dropped, a real-world mistake
    write_dictionary(meta["dictionary"], columns=bad)
    rep = run(meta)
    assert errors_matching(rep, "column mismatch")


def test_invalid_variable_name(meta):
    rows = list(DICT_ROWS)
    rows[1] = ("Country Name", "participant", "radio", "Country",
               "1, Belgium | 2, Chile", "", "y", "", "", "")
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "not a valid REDCap variable name")


def test_invalid_field_type(meta):
    rows = list(DICT_ROWS)
    rows[3] = ("collected_on", "specimen", "datefield", "Collection date",
               "", "date_ymd", "", "", "", "")
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "invalid Field Type")


def test_unknown_validation_type(meta):
    rows = list(DICT_ROWS)
    rows[3] = ("collected_on", "specimen", "text", "Collection date",
               "", "iso_date", "", "", "", "")
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "unknown validation type")


def test_radio_without_choices(meta):
    rows = list(DICT_ROWS)
    rows[1] = ("country", "participant", "radio", "Country", "", "", "y", "", "", "")
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "needs choices")


def test_malformed_choice_syntax(meta):
    rows = list(DICT_ROWS)
    rows[1] = ("country", "participant", "radio", "Country",
               "Belgium | Chile", "", "y", "", "", "")
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "needs choices")


def test_required_flag_must_be_y(meta):
    rows = list(DICT_ROWS)
    rows[0] = ("record_id", "participant", "text", "Sample ID", "", "", "TRUE", "", "", "")
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "Required Field?")


def test_duplicate_variable_name(meta):
    rows = list(DICT_ROWS) + [
        ("country", "specimen", "text", "Country again", "", "", "", "", "", "")
    ]
    write_dictionary(meta["dictionary"], rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "duplicate variable name")


def test_first_records_column_must_be_record_id(meta):
    cols = [RECORD_COLUMNS[1], RECORD_COLUMNS[0]] + RECORD_COLUMNS[2:]
    rows = [[r[1], r[0]] + r[2:] for r in RECORD_ROWS]
    write_records(meta["records"], columns=cols, rows=rows)
    rep = run(meta)
    assert errors_matching(rep, "first column must be the record ID field")


# ------------------------------------------- 5. exit status

def test_main_returns_nonzero_on_error(meta, tmp_path):
    write_records(meta["records"], rows=[
        ["BEHE003", "9", "1", "2019-03-04", "1000", "1", "1", "0"],
        RECORD_ROWS[1],
    ])
    code = vm.main([
        "--dictionary", str(meta["dictionary"]),
        "--records", str(meta["records"]),
        "--biosample", str(meta["biosample"]),
        "--sra", str(meta["sra"]),
        "--manifest", str(meta["manifest"]),
        "--report", str(tmp_path / "report.txt"),
    ])
    assert code == 1
    assert (tmp_path / "report.txt").read_text().strip().endswith("FAIL")


def test_main_returns_zero_when_clean(meta, tmp_path):
    code = vm.main([
        "--dictionary", str(meta["dictionary"]),
        "--records", str(meta["records"]),
        "--biosample", str(meta["biosample"]),
        "--sra", str(meta["sra"]),
        "--manifest", str(meta["manifest"]),
        "--report", str(tmp_path / "report.txt"),
    ])
    assert code == 0
    assert (tmp_path / "report.txt").read_text().strip().endswith("PASS")
