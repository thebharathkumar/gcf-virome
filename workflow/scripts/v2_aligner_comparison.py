"""bowtie2 versus minimap2 on the host-depleted reads.

Outputs:
  per_aligner.tsv             totals per aligner
  per_sample.tsv              viral pairs, references present and wall time
                              per sample and aligner
  breadth_per_reference.tsv   reads and breadth per (sample, reference) for
                              both aligners, every reference either aligner
                              gave at least one filtered read
  agreement.tsv               every (sample, reference) present under at
                              least one aligner, labelled both / bowtie2_only
                              / minimap2_only, plus sample "ANY" rows for the
                              study-wide reference sets
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2lib  # noqa: E402

inp = snakemake.input  # noqa: F821
out = snakemake.output  # noqa: F821

ALIGNERS = ["bowtie2", "minimap2"]


def by_sample(paths, suffix):
    return {os.path.basename(p)[: -len(suffix)]: p for p in paths}


def bench_seconds(path):
    with open(path) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    return float(row["s"]), float(row["max_rss"]) if row.get("max_rss") not in (None, "", "-") else None


with open(inp.samples) as fh:
    samples = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]

with open(inp.contigs) as fh:
    lengths = {r["accession"]: int(r["length"]) for r in csv.DictReader(fh, delimiter="\t")}

flag = {"bowtie2": by_sample(inp.flagstat_bt2, ".flagstat"),
        "minimap2": by_sample(inp.flagstat_mm2, ".flagstat")}
bench = {"bowtie2": by_sample(inp.bench_bt2, ".tsv"),
         "minimap2": by_sample(inp.bench_mm2, ".tsv")}
long = {"bowtie2": v2lib.read_long_table(inp.long_bt2),
        "minimap2": v2lib.read_long_table(inp.long_mm2)}

per_sample = []
for a in ALIGNERS:
    for s in samples:
        with open(flag[a][s]) as fh:
            pairs = v2lib.flagstat_proper_pairs(fh.read())
        secs, rss = bench_seconds(bench[a][s])
        per_sample.append({
            "aligner": a,
            "sample_id": s,
            "viral_pairs": pairs,
            "refs_present": len(v2lib.present_set(long[a], s)),
            "refs_with_any_read": sum(1 for (x, _) in long[a] if x == s),
            "wall_seconds": f"{secs:.2f}",
            "max_rss_mb": f"{rss:.1f}" if rss is not None else "NA",
        })

with open(out.per_sample, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(per_sample[0]), delimiter="\t")
    w.writeheader()
    w.writerows(per_sample)

study = {a: v2lib.present_set(long[a]) for a in ALIGNERS}
agree_study = v2lib.agreement(study["bowtie2"], study["minimap2"])

per_aligner = []
for a in ALIGNERS:
    rows = [r for r in per_sample if r["aligner"] == a]
    per_aligner.append({
        "aligner": a,
        "viral_pairs": sum(r["viral_pairs"] for r in rows),
        "refs_present_any_sample": len(study[a]),
        "sample_ref_calls": sum(r["refs_present"] for r in rows),
        "refs_present_both_aligners": len(agree_study["both"]),
        "refs_present_only_this_aligner": len(
            agree_study["only_a"] if a == "bowtie2" else agree_study["only_b"]
        ),
        "wall_seconds_total": f"{sum(float(r['wall_seconds']) for r in rows):.2f}",
    })

with open(out.per_aligner, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(per_aligner[0]), delimiter="\t")
    w.writeheader()
    w.writerows(per_aligner)

keys = sorted(set(long["bowtie2"]) | set(long["minimap2"]))
with open(out.breadth, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample_id", "accession", "virus", "length",
                "reads_bowtie2", "breadth_bowtie2", "present_bowtie2",
                "reads_minimap2", "breadth_minimap2", "present_minimap2"])
    for k in keys:
        b = long["bowtie2"].get(k)
        m = long["minimap2"].get(k)
        w.writerow([
            k[0], k[1], (b or m)["virus"], lengths.get(k[1], "NA"),
            b["reads"] if b else 0, f"{b['breadth']:.6f}" if b else "0.000000",
            b["present"] if b else 0,
            m["reads"] if m else 0, f"{m['breadth']:.6f}" if m else "0.000000",
            m["present"] if m else 0,
        ])

with open(out.agreement, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample_id", "accession", "virus", "call"])
    for s in samples:
        ag = v2lib.agreement(v2lib.present_set(long["bowtie2"], s),
                             v2lib.present_set(long["minimap2"], s))
        for label, accs in (("both", ag["both"]), ("bowtie2_only", ag["only_a"]),
                            ("minimap2_only", ag["only_b"])):
            for acc in accs:
                r = long["bowtie2"].get((s, acc)) or long["minimap2"].get((s, acc))
                w.writerow([s, acc, r["virus"], label])
    for label, accs in (("both", agree_study["both"]),
                        ("bowtie2_only", agree_study["only_a"]),
                        ("minimap2_only", agree_study["only_b"])):
        for acc in accs:
            r = next(v for (x, y), v in list(long["bowtie2"].items()) + list(long["minimap2"].items())
                     if y == acc)
            w.writerow(["ANY", acc, r["virus"], label])

for r in per_aligner:
    print(r)
