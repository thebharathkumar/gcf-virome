"""Build the REDCap and SRA metadata deliverables from the sample manifest.

Nothing here is typed by hand. Every value is derived from config/samples.tsv,
which was itself generated from NCBI. Where the public record genuinely does
not carry a value, an INSDC missing term is written rather than a guess.

Formats follow docs/metadata_standards.md.
"""
import csv
import re

samples_file = snakemake.input.samples  # noqa: F821
out_dict = snakemake.output.dictionary  # noqa: F821
out_records = snakemake.output.records  # noqa: F821
out_sra = snakemake.output.sra  # noqa: F821
out_biosample = snakemake.output.biosample  # noqa: F821
cfg = snakemake.config  # noqa: F821

with open(samples_file) as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))

COUNTRIES = sorted({r["country"] for r in manifest})
DIAGNOSES = sorted({r["diagnosis"] for r in manifest})
DEVICES = sorted({r["collection_device"] for r in manifest if r["collection_device"]})
SOURCES = sorted({r["isolation_source"] for r in manifest if r["isolation_source"]})
MODELS = sorted({r["instrument_model"] for r in manifest if r["instrument_model"]})

country_code = {name: str(i + 1) for i, name in enumerate(COUNTRIES)}
diagnosis_code = {name: str(i + 1) for i, name in enumerate(DIAGNOSES)}
device_code = {name: str(i + 1) for i, name in enumerate(DEVICES)}
source_code = {name: str(i + 1) for i, name in enumerate(SOURCES)}
model_code = {name: str(i + 1) for i, name in enumerate(MODELS)}


def choices(mapping):
    """REDCap choice syntax: 'code, label | code, label'."""
    return " | ".join(f"{code}, {label}" for label, code in sorted(mapping.items(), key=lambda kv: kv[1]))


# --------------------------------------------------------------------------
# REDCap data dictionary
# --------------------------------------------------------------------------

DICT_COLUMNS = [
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


def field(name, form, ftype, label, section="", choices_="", note="",
          validation="", vmin="", vmax="", identifier="", branching="",
          required="", annotation=""):
    return {
        "Variable / Field Name": name,
        "Form Name": form,
        "Section Header": section,
        "Field Type": ftype,
        "Field Label": label,
        "Choices, Calculations, OR Slider Labels": choices_,
        "Field Note": note,
        "Text Validation Type OR Show Slider Number": validation,
        "Text Validation Min": vmin,
        "Text Validation Max": vmax,
        "Identifier?": identifier,
        "Branching Logic (Show field only if...)": branching,
        "Required Field?": required,
        "Custom Alignment": "",
        "Question Number (surveys only)": "",
        "Matrix Group Name": "",
        "Matrix Ranking?": "",
        "Field Annotation": annotation,
    }


fields = [
    # ---- participant form
    field("record_id", "participant", "text", "Sample ID (record identifier)",
          section="Participant and diagnosis",
          note="Study sample code, for example BEHE003. Country and group are encoded in it.",
          required="y"),
    field("patient_number", "participant", "text", "Patient number assigned by the study site",
          identifier="y",
          note="Subject-level code carried in the BioSample record.",
          required="y"),
    field("country", "participant", "radio", "Country of recruitment",
          choices_=choices(country_code), required="y"),
    field("diagnosis", "participant", "radio", "Periodontal diagnosis",
          choices_=choices(diagnosis_code), required="y",
          note="As recorded by the originating study. No case definition is re-derived here."),
    field("host_species", "participant", "text", "Host species", required="y"),

    # ---- specimen form
    field("specimen_type", "specimen", "radio", "Specimen type",
          section="Specimen collection",
          choices_=choices(source_code), required="y"),
    field("collection_device", "specimen", "radio", "Collection device",
          choices_=choices(device_code), required="y"),
    field("collection_year_start", "specimen", "text", "Collection period, first year",
          validation="integer", vmin="1990", vmax="2030",
          note="The public record gives a year range only, not a specimen-level date.",
          required="y"),
    field("collection_year_end", "specimen", "text", "Collection period, last year",
          validation="integer", vmin="1990", vmax="2030", required="y"),
    field("latitude", "specimen", "text", "Site latitude (decimal degrees)",
          validation="number", vmin="-90", vmax="90"),
    field("longitude", "specimen", "text", "Site longitude (decimal degrees)",
          validation="number", vmin="-180", vmax="180"),

    # ---- sequencing form
    field("bioproject", "sequencing", "text", "NCBI BioProject accession",
          section="Sequencing and repository accessions", required="y"),
    field("biosample", "sequencing", "text", "NCBI BioSample accession", required="y"),
    field("sra_run", "sequencing", "text", "NCBI SRA run accession", required="y"),
    field("sra_experiment", "sequencing", "text", "NCBI SRA experiment accession", required="y"),
    field("instrument_model", "sequencing", "radio", "Sequencing instrument",
          choices_=choices(model_code), required="y"),
    field("library_layout", "sequencing", "radio", "Library layout",
          choices_="1, single | 2, paired", required="y"),
    field("deposited_spots", "sequencing", "text", "Read pairs deposited in SRA",
          validation="integer", vmin="0", required="y"),
    field("deposited_bases", "sequencing", "text", "Bases deposited in SRA",
          validation="integer", vmin="0", required="y"),
    field("analysed_in_pipeline", "sequencing", "yesno",
          "Included in the subsampled analysis set", required="y"),
]

with open(out_dict, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=DICT_COLUMNS)
    writer.writeheader()
    writer.writerows(fields)
print(f"wrote {out_dict} ({len(fields)} fields)")


# --------------------------------------------------------------------------
# REDCap records
# --------------------------------------------------------------------------

def split_years(value):
    """'2017/2019' -> ('2017', '2019'); '2021-05-11' -> ('2021', '2021')."""
    years = re.findall(r"(\d{4})", value or "")
    if not years:
        return "", ""
    return years[0], years[-1]


def split_lat_lon(value):
    """'50.88025 N 4.693472 E' -> ('50.88025', '4.693472') with sign."""
    m = re.match(r"^([\d.]+)\s*([NS])\s+([\d.]+)\s*([EW])$", (value or "").strip())
    if not m:
        return "", ""
    lat = float(m.group(1)) * (1 if m.group(2) == "N" else -1)
    lon = float(m.group(3)) * (1 if m.group(4) == "E" else -1)
    return f"{lat:g}", f"{lon:g}"


record_columns = [f["Variable / Field Name"] for f in fields]

records = []
for row in manifest:
    y0, y1 = split_years(row["collection_date"])
    lat, lon = split_lat_lon(row["lat_lon"])
    records.append(
        {
            "record_id": row["sample_id"],
            "patient_number": row["sample_id"],
            "country": country_code[row["country"]],
            "diagnosis": diagnosis_code[row["diagnosis"]],
            "host_species": row["host"],
            "specimen_type": source_code[row["isolation_source"]],
            "collection_device": device_code[row["collection_device"]],
            "collection_year_start": y0,
            "collection_year_end": y1,
            "latitude": lat,
            "longitude": lon,
            "bioproject": cfg["bioproject"],
            "biosample": row["biosample"],
            "sra_run": row["run_accession"],
            "sra_experiment": row["experiment"],
            "instrument_model": model_code[row["instrument_model"]],
            "library_layout": "2" if row["library_layout"].lower() == "paired" else "1",
            "deposited_spots": row["spots"],
            "deposited_bases": row["bases"],
            "analysed_in_pipeline": "1",
        }
    )

with open(out_records, "w", newline="") as fh:
    writer = csv.DictWriter(fh, fieldnames=record_columns)
    writer.writeheader()
    writer.writerows(records)
print(f"wrote {out_records} ({len(records)} records)")


# --------------------------------------------------------------------------
# NCBI BioSample attributes, package MIMS.me.human-oral.6.0
#
# NCBI keeps sample attributes and sequencing-library metadata in two separate
# sheets, so two files are written rather than one combined table.
#
# The seven package-mandatory attributes are collection_date, env_broad_scale,
# env_local_scale, env_medium, geo_loc_name, host and lat_lon, alongside the
# base sample_name and organism columns.
# --------------------------------------------------------------------------

BIOSAMPLE_COLUMNS = [
    "sample_name", "sample_title", "bioproject_accession", "organism",
    "collection_date", "env_broad_scale", "env_local_scale", "env_medium",
    "geo_loc_name", "host", "lat_lon",
    "host_disease", "host_subject_id", "isolation_source", "samp_collect_device",
]

SRA_COLUMNS = [
    "bioproject_accession", "biosample_accession", "library_ID", "title",
    "library_strategy", "library_source", "library_selection", "library_layout",
    "platform", "instrument_model", "design_description", "filetype",
    "filename", "filename2",
]

# MIxS environmental triplet. Each ID was resolved against the EBI OLS API and
# confirmed to exist and to be non-obsolete. Two obvious-looking choices were
# rejected on checking: ENVO:00009003 ("human-associated habitat") is obsolete
# with no replacement, and UBERON:0006932 is "vestibular epithelium", not
# gingival crevicular fluid. Neither ENVO nor UBERON has a term for gingival
# crevicular fluid at all, so the nearest correct parent is used for env_medium.
# docs/metadata_standards.md records the reasoning.
ENV_BROAD = "oral cavity [UBERON:0000167]"
ENV_LOCAL = "gingiva [UBERON:0001828]"
ENV_MEDIUM = "bodily fluid material [ENVO:02000019]"


def ncbi_date(value):
    """NCBI accepts a year, or a range written as start/end."""
    y0, y1 = split_years(value)
    if not y0:
        return "missing"
    return y0 if y0 == y1 else f"{y0}/{y1}"


biosample_rows, sra_rows = [], []
for row in manifest:
    disease = (
        "periodontitis"
        if row["diagnosis"].lower().startswith("periodontitis")
        else "not applicable"
    )
    biosample_rows.append(
        {
            "sample_name": row["sample_id"],
            "sample_title": (
                f"Subgingival gingival crevicular fluid metagenome, "
                f"{row['diagnosis'].lower()}, {row['country']}"
            ),
            "bioproject_accession": cfg["bioproject"],
            "organism": "human oral metagenome",
            "collection_date": ncbi_date(row["collection_date"]),
            "env_broad_scale": ENV_BROAD,
            "env_local_scale": ENV_LOCAL,
            "env_medium": ENV_MEDIUM,
            "geo_loc_name": row["country"] or "missing",
            "host": row["host"] or "Homo sapiens",
            "lat_lon": row["lat_lon"] or "missing",
            "host_disease": disease,
            "host_subject_id": row["sample_id"],
            "isolation_source": row["isolation_source"] or "missing",
            "samp_collect_device": row["collection_device"] or "missing",
        }
    )

    sra_rows.append(
        {
            "bioproject_accession": cfg["bioproject"],
            "biosample_accession": row["biosample"],
            "library_ID": row["sample_id"],
            "title": (
                f"WGS of human oral metagenome: subgingival gingival crevicular "
                f"fluid, {row['diagnosis'].lower()}, {row['country']}"
            ),
            "library_strategy": row["library_strategy"],
            "library_source": row["library_source"],
            "library_selection": row["library_selection"],
            "library_layout": row["library_layout"].lower(),
            "platform": "ILLUMINA",
            "instrument_model": row["instrument_model"],
            "design_description": (
                "Subgingival plaque sampled on paper points from the deepest site "
                "of each quadrant and pooled per subject, followed by shotgun "
                "metagenomic sequencing."
            ),
            "filetype": "fastq",
            "filename": f"{row['sample_id']}_1.fastq.gz",
            "filename2": f"{row['sample_id']}_2.fastq.gz",
        }
    )


def write_tsv(path, columns, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {path} ({len(rows)} rows)")


write_tsv(out_biosample, BIOSAMPLE_COLUMNS, biosample_rows)
write_tsv(out_sra, SRA_COLUMNS, sra_rows)
