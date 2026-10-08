"""Explain aligner disagreement from the reads themselves.

For each (sample, reference) present under one aligner and not the other:

  1. take every read the detecting aligner counted on that reference (filtered
     BAM, primary records), and
  2. look the same read and mate up in the other aligner's unfiltered BAM,
     then classify what the other aligner did with it (see
     v2lib.classify_in_other).

discordant_reads.tsv has one row per read. discordant_summary.tsv has one row
per discordant (sample, reference) with the category counts and the most
common alternative reference, which is what the README discussion is built on.
"""
import csv
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2lib  # noqa: E402

inp = snakemake.input  # noqa: F821
out = snakemake.output  # noqa: F821
base = snakemake.params.dir  # noqa: F821
min_mapq = int(snakemake.params.mapq)  # noqa: F821

RUN = {"bowtie2": "dehost_bowtie2", "minimap2": "dehost_minimap2"}
OTHER = {"bowtie2": "minimap2", "minimap2": "bowtie2"}

with open(inp.contigs) as fh:
    desc = {r["accession"]: r["description"] for r in csv.DictReader(fh, delimiter="\t")}


def sam(args):
    res = subprocess.run(["samtools", "view"] + args, check=True,
                         capture_output=True, text=True)
    return [v2lib.parse_sam_line(x) for x in res.stdout.splitlines() if x]


discordant = []
with open(inp.agreement) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["sample_id"] == "ANY" or r["call"] == "both":
            continue
        discordant.append((r["sample_id"], r["accession"], r["call"].split("_")[0]))

read_rows = []
summary_rows = []
for sample, acc, found_by in discordant:
    other = OTHER[found_by]
    filt = f"{base}/{RUN[found_by]}/bam/{sample}.bam"
    raw_other = f"{base}/{RUN[other]}/raw/{sample}.bam"
    raw_self = f"{base}/{RUN[found_by]}/raw/{sample}.bam"

    counted = [r for r in sam(["-F", "0x900", filt, acc])]
    names = sorted({r["qname"] for r in counted})
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
        tf.write("\n".join(names) + "\n")
        name_file = tf.name
    try:
        other_recs = defaultdict(list)
        for rec in sam(["-N", name_file, raw_other]):
            other_recs[(rec["qname"], rec["mate"])].append(rec)
        # Reads on this reference in the other aligner's output that it did
        # not count, so the summary can show what fell just short.
        other_on_ref = sam([raw_other, acc])
    finally:
        os.unlink(name_file)

    cats = Counter()
    alt_refs = Counter()
    for rec in counted:
        o = other_recs.get((rec["qname"], rec["mate"]), [])
        cat = v2lib.classify_in_other(acc, o, min_mapq)
        cats[cat] += 1
        prim = next((x for x in o if not x["flag"] & 0x904), None)
        if cat == "other_ref" and prim:
            alt_refs[prim["rname"]] += 1
        read_rows.append({
            "sample_id": sample,
            "accession": acc,
            "found_by": found_by,
            "qname": rec["qname"],
            "mate": rec["mate"],
            "pos": rec["pos"],
            "mapq": rec["mapq"],
            "cigar": rec["cigar"],
            "nm": rec["nm"] if rec["nm"] is not None else "NA",
            "other_category": cat,
            "other_rname": prim["rname"] if prim else "*",
            "other_pos": prim["pos"] if prim else 0,
            "other_mapq": prim["mapq"] if prim else "NA",
            "other_cigar": prim["cigar"] if prim else "*",
            "other_nm": prim["nm"] if prim and prim["nm"] is not None else "NA",
            "other_proper_pair": int(bool(prim and prim["flag"] & 0x2)),
        })

    other_primary = [x for x in other_on_ref if not x["flag"] & 0x904]
    other_pass = [x for x in other_primary if v2lib.passes_filters(x, min_mapq)]
    mq = [r["mapq"] for r in counted]
    nm = [r["nm"] for r in counted if r["nm"] is not None]
    top_alt = alt_refs.most_common(1)
    summary_rows.append({
        "sample_id": sample,
        "accession": acc,
        "virus": desc.get(acc, acc),
        "found_by": found_by,
        "reads_counted": len(counted),
        "median_mapq": sorted(mq)[len(mq) // 2] if mq else "NA",
        "mean_nm": f"{sum(nm) / len(nm):.2f}" if nm else "NA",
        "other_unaligned": cats["unaligned"],
        "other_other_ref": cats["other_ref"],
        "other_same_ref_low_mapq": cats["same_ref_low_mapq"],
        "other_same_ref_not_proper": cats["same_ref_not_proper"],
        "other_same_ref_pass": cats["same_ref_pass"],
        "top_alt_ref": top_alt[0][0] if top_alt else "",
        "top_alt_ref_reads": top_alt[0][1] if top_alt else 0,
        "top_alt_ref_desc": desc.get(top_alt[0][0], "") if top_alt else "",
        "other_primary_on_ref_any_mapq": len(other_primary),
        "other_passing_on_ref": len(other_pass),
        "self_raw_primary_on_ref": len(
            [x for x in sam([raw_self, acc]) if not x["flag"] & 0x904]
        ),
    })

read_cols = ["sample_id", "accession", "found_by", "qname", "mate", "pos", "mapq",
             "cigar", "nm", "other_category", "other_rname", "other_pos",
             "other_mapq", "other_cigar", "other_nm", "other_proper_pair"]
with open(out.reads, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=read_cols, delimiter="\t")
    w.writeheader()
    w.writerows(read_rows)

sum_cols = ["sample_id", "accession", "virus", "found_by", "reads_counted",
            "median_mapq", "mean_nm", "other_unaligned", "other_other_ref",
            "other_same_ref_low_mapq", "other_same_ref_not_proper",
            "other_same_ref_pass", "top_alt_ref", "top_alt_ref_reads",
            "top_alt_ref_desc", "other_primary_on_ref_any_mapq",
            "other_passing_on_ref", "self_raw_primary_on_ref"]
with open(out.summary, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=sum_cols, delimiter="\t")
    w.writeheader()
    w.writerows(summary_rows)

print(f"{len(discordant)} discordant (sample, reference) calls, "
      f"{len(read_rows)} reads examined")
