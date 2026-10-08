"""Per-sample effect of host read removal on the viral alignment.

Writes two tables:

  host_depletion_per_sample.tsv  one row per sample plus a TOTAL row
  reference_changes.tsv          every (sample, reference) whose presence call
                                 changed between the baseline and the
                                 host-depleted alignment
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2lib  # noqa: E402

inp = snakemake.input  # noqa: F821
out = snakemake.output  # noqa: F821
herv = snakemake.params.herv  # noqa: F821


def by_sample(paths, suffix):
    return {os.path.basename(p)[: -len(suffix)]: p for p in paths}


def text(path):
    with open(path) as fh:
        return fh.read()


with open(inp.samples) as fh:
    samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]

fastp = by_sample(inp.fastp, ".json")
hostile = by_sample(inp.hostile, ".json")
fs_before = by_sample(inp.flagstat_before, ".flagstat")
fs_after = by_sample(inp.flagstat_after, ".flagstat")
before = v2lib.read_long_table(inp.long_before)
after = v2lib.read_long_table(inp.long_after)

rows = []
for s in samples:
    rows.append(
        v2lib.host_depletion_row(
            sample=s,
            post_qc_pairs=v2lib.fastp_post_qc_pairs(text(fastp[s])),
            hostile=v2lib.hostile_stats(text(hostile[s])),
            before=before,
            after=after,
            accession=herv,
            pairs_before=v2lib.flagstat_proper_pairs(text(fs_before[s])),
            pairs_after=v2lib.flagstat_proper_pairs(text(fs_after[s])),
        )
    )

total = {"sample_id": "TOTAL"}
for k in rows[0]:
    if k == "sample_id":
        continue
    if k == "host_pct_of_post_qc":
        continue
    total[k] = sum(r[k] for r in rows)
total["input_matches_post_qc"] = int(all(r["input_matches_post_qc"] for r in rows))
total["host_pct_of_post_qc"] = round(
    100.0 * total["host_reads_removed"] / (2 * total["post_qc_pairs"]), 4
)
# Distinct references present in any sample, not a sum over samples.
total["refs_present_before"] = len(v2lib.present_set(before))
total["refs_present_after"] = len(v2lib.present_set(after))

fields = list(rows[0])
with open(out.per_sample, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(rows)
    w.writerow({k: total[k] for k in fields})

changes = []
for key in sorted(set(before) | set(after)):
    b, a = before.get(key), after.get(key)
    pb = b["present"] if b else 0
    pa = a["present"] if a else 0
    rb = b["reads"] if b else 0
    ra = a["reads"] if a else 0
    if pb != pa or rb != ra:
        changes.append({
            "sample_id": key[0],
            "accession": key[1],
            "virus": (b or a)["virus"],
            "reads_before": rb,
            "reads_after": ra,
            "breadth_before": f"{b['breadth']:.6f}" if b else "0.000000",
            "breadth_after": f"{a['breadth']:.6f}" if a else "0.000000",
            "present_before": pb,
            "present_after": pa,
        })

with open(out.reference_changes, "w", newline="") as fh:
    cols = ["sample_id", "accession", "virus", "reads_before", "reads_after",
            "breadth_before", "breadth_after", "present_before", "present_after"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
    w.writeheader()
    w.writerows(changes)

print(f"{len(rows)} samples; host reads removed {total['host_reads_removed']:,} "
      f"({total['host_pct_of_post_qc']}% of post-QC reads)")
print(f"viral pairs {total['viral_pairs_before']:,} -> {total['viral_pairs_after']:,}")
print(f"HERV-K113 reads {total['herv_k113_reads_before']:,} -> "
      f"{total['herv_k113_reads_after']:,}")
print(f"{len(changes)} (sample, reference) rows changed")
