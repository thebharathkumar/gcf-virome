"""Check that the v2 baseline alignment reproduces the committed v1 run.

The "before" side of the host depletion experiment is realigned inside the v2
tree so that before and after share one software build. This compares it with
the v1 tables committed to the repository. The output lists every difference;
a single PASS row means there were none.
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import v2lib  # noqa: E402

inp = snakemake.input  # noqa: F821
out = snakemake.output[0]  # noqa: F821

v1 = v2lib.read_long_table(inp.v1_long)
v2 = v2lib.read_long_table(inp.v2_long)
diffs = v2lib.compare_long_tables(v1, v2)

pair_diffs = []
with open(inp.v1_retention) as fh:
    v1_pairs = {r["sample_id"]: int(r["mapped_pairs"]) for r in csv.DictReader(fh, delimiter="\t")}
for p in inp.v2_flagstat:
    s = os.path.basename(p)[: -len(".flagstat")]
    with open(p) as fh:
        n = v2lib.flagstat_proper_pairs(fh.read())
    if v1_pairs.get(s) != n:
        pair_diffs.append((s, v1_pairs.get(s), n))

cols = ["check", "sample_id", "accession", "v1", "v2"]
with open(out, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(cols)
    for d in diffs:
        w.writerow(["reads", d["sample_id"], d["accession"], d["reads_v1"], d["reads_v2"]])
        w.writerow(["breadth", d["sample_id"], d["accession"], d["breadth_v1"], d["breadth_v2"]])
        w.writerow(["present", d["sample_id"], d["accession"], d["present_v1"], d["present_v2"]])
    for s, a, b in pair_diffs:
        w.writerow(["viral_pairs", s, "", a, b])
    if not diffs and not pair_diffs:
        w.writerow(["PASS", "all", "all", "identical", "identical"])

print(f"{len(diffs)} count-table differences, {len(pair_diffs)} viral-pair differences")
