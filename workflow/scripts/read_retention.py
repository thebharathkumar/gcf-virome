"""Read counts retained at every step, one row per sample.

Columns:
  sra_spots         read pairs deposited in SRA, from the NCBI run table
  subsampled_pairs  read pairs after seqtk, the pipeline's actual input
  post_qc_pairs     read pairs surviving fastp
  qc_retained_pct   post_qc_pairs / subsampled_pairs
  mapped_pairs      properly paired alignments to RefSeq viral above the MAPQ
                    cutoff, divided by two to express as pairs
  mapped_pct        mapped_pairs / post_qc_pairs
"""
import csv
import json
import os
import re

fastp_files = snakemake.input.fastp  # noqa: F821
flagstat_files = snakemake.input.flagstat  # noqa: F821
samples_file = snakemake.input.samples  # noqa: F821
out = snakemake.output[0]  # noqa: F821

with open(samples_file) as fh:
    manifest = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}

fastp_by_sample = {os.path.basename(p)[: -len(".json")]: p for p in fastp_files}
flagstat_by_sample = {os.path.basename(p)[: -len(".flagstat")]: p for p in flagstat_files}


def read_fastp(path):
    with open(path) as fh:
        d = json.load(fh)
    before = d["summary"]["before_filtering"]["total_reads"]
    after = d["summary"]["after_filtering"]["total_reads"]
    # fastp counts reads, not pairs, for paired input.
    return before // 2, after // 2


def read_flagstat(path):
    """Pull 'properly paired' read count from samtools flagstat."""
    with open(path) as fh:
        text = fh.read()
    m = re.search(r"^(\d+) \+ \d+ properly paired", text, re.M)
    if m:
        return int(m.group(1)) // 2
    m = re.search(r"^(\d+) \+ \d+ in total", text, re.M)
    return (int(m.group(1)) // 2) if m else 0


rows = []
for sample_id, rec in manifest.items():
    sub_pairs, post_pairs = read_fastp(fastp_by_sample[sample_id])
    mapped_pairs = read_flagstat(flagstat_by_sample[sample_id])
    rows.append(
        {
            "sample_id": sample_id,
            "run_accession": rec["run_accession"],
            "diagnosis": rec["diagnosis"],
            "country": rec["country"],
            "sra_spots": int(rec["spots"]),
            "subsampled_pairs": sub_pairs,
            "post_qc_pairs": post_pairs,
            "qc_retained_pct": round(100.0 * post_pairs / sub_pairs, 2) if sub_pairs else 0.0,
            "mapped_pairs": mapped_pairs,
            "mapped_pct": round(100.0 * mapped_pairs / post_pairs, 4) if post_pairs else 0.0,
        }
    )

rows.sort(key=lambda r: r["sample_id"])

with open(out, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader()
    w.writerows(rows)

total_sub = sum(r["subsampled_pairs"] for r in rows)
total_post = sum(r["post_qc_pairs"] for r in rows)
total_map = sum(r["mapped_pairs"] for r in rows)
print(f"{len(rows)} samples")
print(f"subsampled pairs : {total_sub:,}")
print(f"post-QC pairs    : {total_post:,} ({100.0 * total_post / total_sub:.2f}%)")
print(f"viral mapped     : {total_map:,} ({100.0 * total_map / total_post:.4f}% of post-QC)")
