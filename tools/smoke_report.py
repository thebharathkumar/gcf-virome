"""Exercise build_count_matrix.py and make_report.py on synthetic inputs.

Values are synthetic and land in build/, never results/. This exists so that a
formatting or layout bug surfaces before the real run finishes.
"""
import builtins
import csv
import random
import runpy
import types
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build" / "smoke_report"
DIV = ROOT / "build" / "smoke_diversity"
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "coverage").mkdir(exist_ok=True)

random.seed(11)
cfg = yaml.safe_load((ROOT / "config" / "config.yaml").read_text())

with open(ROOT / "config" / "samples.tsv") as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))
samples = [r["sample_id"] for r in manifest]

# ---- fake reference contig table
contigs = OUT / "contigs.tsv"
refs = [(f"NC_{900000 + i}.1", i) for i in range(40)]
with open(contigs, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["accession", "length", "description"])
    for acc, i in refs:
        w.writerow([acc, 40000 + i * 100, f"Fake virus {i} isolate X, complete genome"])

# ---- fake samtools coverage output per sample
for s in samples:
    p = OUT / "coverage" / f"{s}.coverage.tsv"
    with open(p, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["#rname", "startpos", "endpos", "numreads", "covbases",
                    "coverage", "meandepth", "meanbaseq", "meanmapq"])
        for acc, i in refs:
            if random.random() < 0.4:
                continue
            nreads = random.randint(1, 500)
            cov = round(random.uniform(0.05, 40.0), 4)
            w.writerow([acc, 1, 40000 + i * 100, nreads, 100, cov, 1.2, 36, 42])

fake = types.SimpleNamespace(
    input=types.SimpleNamespace(
        coverage=[str(OUT / "coverage" / f"{s}.coverage.tsv") for s in samples],
        contigs=str(contigs),
        samples=str(ROOT / "config" / "samples.tsv"),
    ),
    output=types.SimpleNamespace(
        counts=str(OUT / "virus_count_matrix.tsv"),
        presence=str(OUT / "virus_presence_matrix.tsv"),
        long=str(OUT / "virus_counts_long.tsv"),
    ),
    params=types.SimpleNamespace(
        min_reads=cfg["presence"]["min_reads"],
        min_breadth=cfg["presence"]["min_breadth"],
    ),
    config=cfg,
)
builtins.snakemake = fake
runpy.run_path(str(ROOT / "workflow" / "scripts" / "build_count_matrix.py"), run_name="__main__")

# ---- fake read retention table
retention = OUT / "read_retention.tsv"
with open(retention, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample_id", "run_accession", "diagnosis", "country", "sra_spots",
                "subsampled_pairs", "post_qc_pairs", "qc_retained_pct",
                "mapped_pairs", "mapped_pct"])
    for r in manifest:
        sub = cfg["subsample"]["n_pairs"]
        post = int(sub * random.uniform(0.93, 0.99))
        mapped = random.randint(20, 900)
        w.writerow([r["sample_id"], r["run_accession"], r["diagnosis"], r["country"],
                    r["spots"], sub, post, round(100 * post / sub, 2),
                    mapped, round(100 * mapped / post, 4)])

# ---- fake reference provenance
refmeta = OUT / "refseq_viral.source.txt"
refmeta.write_text(
    "url\thttps://example.invalid/viral.fna.gz\n"
    "refseq_release\t237\nsha256\tdeadbeef\nbytes\t178194466\n"
    "retrieved\t2026-09-27T00:00:00Z\n"
)

fake2 = types.SimpleNamespace(
    input=types.SimpleNamespace(
        retention=str(retention),
        alpha=str(DIV / "alpha_diversity.tsv"),
        alpha_tests=str(DIV / "alpha_group_tests.tsv"),
        permanova=str(DIV / "permanova.tsv"),
        pcoa=str(DIV / "pcoa_coordinates.tsv"),
        counts=str(OUT / "virus_count_matrix.tsv"),
        presence=str(OUT / "virus_presence_matrix.tsv"),
        figure=str(DIV / "diversity.png"),
        reference=str(refmeta),
        samples=str(ROOT / "config" / "samples.tsv"),
    ),
    output=[str(OUT / "report.pdf")],
    config=cfg,
)
builtins.snakemake = fake2
runpy.run_path(str(ROOT / "workflow" / "scripts" / "make_report.py"), run_name="__main__")

pdf = OUT / "report.pdf"
print(f"\nreport.pdf exists: {pdf.exists()}  bytes: {pdf.stat().st_size if pdf.exists() else 0}")
