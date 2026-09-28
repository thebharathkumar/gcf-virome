"""Exercise R/diversity.R on a synthetic count matrix.

This only checks that the code path runs and writes every expected file. The
numbers it produces are meaningless and are written to build/, never to
results/. It exists so that an R API problem surfaces now rather than at the
end of a three-hour run.
"""
import csv
import random
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build" / "smoke_diversity"
OUT.mkdir(parents=True, exist_ok=True)

random.seed(7)

with open(ROOT / "config" / "samples.tsv") as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))
samples = [r["sample_id"] for r in manifest]

# 30 fake viruses, a couple of them segmented so the aggregation-by-name path
# is actually taken.
viruses = []
for i in range(30):
    name = f"Fake phage {i // 2}" if i < 6 else f"Fake virus {i}"
    viruses.append((f"NC_{900000 + i}.1", name))

counts_path = OUT / "counts.tsv"
presence_path = OUT / "presence.tsv"

with open(counts_path, "w", newline="") as cf, open(presence_path, "w", newline="") as pf:
    cw = csv.writer(cf, delimiter="\t")
    pw = csv.writer(pf, delimiter="\t")
    cw.writerow(["accession", "virus"] + samples)
    pw.writerow(["accession", "virus"] + samples)
    for acc, name in viruses:
        counts, presence = [], []
        for s in samples:
            if random.random() < 0.45:
                n = random.randint(0, 400)
            else:
                n = 0
            counts.append(n)
            presence.append(1 if n >= 10 else 0)
        cw.writerow([acc, name] + counts)
        pw.writerow([acc, name] + presence)

args = [
    "Rscript", str(ROOT / "R" / "diversity.R"),
    str(counts_path), str(presence_path), str(ROOT / "config" / "samples.tsv"),
    "diagnosis",
    str(OUT / "alpha_diversity.tsv"),
    str(OUT / "alpha_group_tests.tsv"),
    str(OUT / "permanova.tsv"),
    str(OUT / "pcoa_coordinates.tsv"),
    str(OUT / "diversity.png"),
    str(OUT / "r_session_info.txt"),
]

proc = subprocess.run(args, capture_output=True, text=True, cwd=ROOT)
print(proc.stdout)
if proc.returncode != 0:
    print("--- STDERR ---", file=sys.stderr)
    print(proc.stderr, file=sys.stderr)
    sys.exit(proc.returncode)

missing = [
    p.name for p in [
        OUT / "alpha_diversity.tsv", OUT / "alpha_group_tests.tsv",
        OUT / "permanova.tsv", OUT / "pcoa_coordinates.tsv",
        OUT / "diversity.png", OUT / "r_session_info.txt",
    ] if not p.exists()
]
if missing:
    print(f"MISSING OUTPUTS: {missing}", file=sys.stderr)
    sys.exit(1)

print("all expected outputs written (values are synthetic and meaningless)")
if proc.stderr.strip():
    print("--- R stderr (warnings) ---")
    print(proc.stderr[:2000])
