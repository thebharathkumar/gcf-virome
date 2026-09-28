"""Build config/samples.tsv for PRJNA1183294 from live NCBI metadata.

Selection rule, fixed in advance so it cannot be tuned to a result:
  stratify by (country, diagnosis), sort sample names lexicographically,
  take the first N_PER_STRATUM from each stratum.

The full 80-sample table is written alongside it so the selection is auditable.
Run this once; the committed TSV is what the workflow reads.
"""
import collections
import csv
import io
import re
import sys
from pathlib import Path

sys.path.insert(0, "tools")
import ncbi  # noqa: E402

BIOPROJECT = "PRJNA1183294"
N_PER_STRATUM = 2
OUT_DIR = Path("config")
SELECTED = OUT_DIR / "samples.tsv"
FULL = OUT_DIR / "bioproject_all_samples.tsv"

FIELDS = [
    "sample_id",
    "run_accession",
    "biosample",
    "experiment",
    "diagnosis",
    "country",
    "isolation_source",
    "collection_device",
    "collection_date",
    "lat_lon",
    "host",
    "library_strategy",
    "library_source",
    "library_selection",
    "library_layout",
    "instrument_model",
    "spots",
    "bases",
    "size_mb",
]


def fetch_runs():
    _, text = ncbi.runinfo(f"{BIOPROJECT}[BioProject]")
    rows = [r for r in csv.DictReader(io.StringIO(text)) if r.get("Run")]
    if not rows:
        raise SystemExit("no SRA runs returned; refusing to write a manifest")
    return rows


def fetch_biosamples(accessions):
    out = {}
    for i in range(0, len(accessions), 20):
        xml = ncbi.efetch("biosample", accessions[i : i + 20], "full", "xml")
        for block in re.findall(r"<BioSample .*?</BioSample>", xml, re.S):
            acc = re.search(r'accession="([^"]+)"', block)
            if not acc:
                continue
            out[acc.group(1)] = {
                k: re.sub(r"\s+", " ", v).strip()
                for k, v in re.findall(
                    r'<Attribute attribute_name="([^"]+)"[^>]*>(.*?)</Attribute>',
                    block,
                    re.S,
                )
            }
    return out


def build_records(runs, samples):
    records = []
    for r in runs:
        attr = samples.get(r["BioSample"], {})
        diagnosis = attr.get("Diagnosis", "")
        if not diagnosis:
            raise SystemExit(f"{r['BioSample']} has no Diagnosis attribute; aborting")
        records.append(
            {
                "sample_id": r["SampleName"],
                "run_accession": r["Run"],
                "biosample": r["BioSample"],
                "experiment": r["Experiment"],
                "diagnosis": diagnosis,
                "country": attr.get("geo_loc_name", ""),
                "isolation_source": attr.get("isolation_source", ""),
                "collection_device": attr.get("samp_collect_device", ""),
                "collection_date": attr.get("collection_date", ""),
                "lat_lon": attr.get("lat_lon", ""),
                "host": attr.get("host", ""),
                "library_strategy": r["LibraryStrategy"],
                "library_source": r["LibrarySource"],
                "library_selection": r["LibrarySelection"],
                "library_layout": r["LibraryLayout"],
                "instrument_model": r["Model"],
                "spots": r["spots"],
                "bases": r["bases"],
                "size_mb": r["size_MB"],
            }
        )
    return sorted(records, key=lambda d: d["sample_id"])


def select(records):
    strata = collections.defaultdict(list)
    for rec in records:
        strata[(rec["country"], rec["diagnosis"])].append(rec)
    chosen = []
    for key in sorted(strata):
        group = sorted(strata[key], key=lambda d: d["sample_id"])
        if len(group) < N_PER_STRATUM:
            raise SystemExit(f"stratum {key} has only {len(group)} samples")
        chosen.extend(group[:N_PER_STRATUM])
    return sorted(chosen, key=lambda d: d["sample_id"])


def write(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(records)
    print(f"wrote {path} ({len(records)} rows)")


runs = fetch_runs()
samples = fetch_biosamples(sorted({r["BioSample"] for r in runs}))
records = build_records(runs, samples)
write(FULL, records)

chosen = select(records)
write(SELECTED, chosen)

print("\nselected samples:")
total_mb = 0
for rec in chosen:
    total_mb += int(rec["size_mb"])
    print(
        f"  {rec['sample_id']:<9} {rec['run_accession']}  "
        f"{rec['diagnosis']:<22} {rec['country']:<8} "
        f"{int(rec['spots']):>11,} spots  {rec['size_mb']:>5} MB"
    )
print(f"\ntotal download: {total_mb:,} MB ({total_mb / 1024:.1f} GB)")
print("strata:", dict(collections.Counter((r["country"], r["diagnosis"]) for r in chosen)))
