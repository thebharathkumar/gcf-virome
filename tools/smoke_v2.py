"""Run the results/v2 rules end to end on a tiny synthetic dataset.

Checks that every v2 rule runs with the real tools (Hostile, bowtie2,
minimap2, samtools, R) and writes every expected file, and that the summary
scripts agree with what was planted. Nothing here touches results/: the run
happens in a sandbox under build/smoke_v2/ that links to the real workflow
code. The numbers it produces describe synthetic reads and mean nothing.

What is planted, so the outputs can be checked against known answers:

  * a fake "human" genome that contains a copy of a fake HERV sequence, used
    as the Hostile index, so host depletion must remove the HERV reads;
  * a fake viral reference set holding that HERV sequence plus three viruses;
  * per sample, read pairs drawn from the host, from the viruses, and random.

Usage (inside the oral-virome environment):
    python tools/smoke_v2.py            # uses --use-conda for Hostile
    python tools/smoke_v2.py --no-conda # Hostile already on PATH
"""
import argparse
import csv
import gzip
import json
import os
import random
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SB = ROOT / "build" / "smoke_v2"

ap = argparse.ArgumentParser()
ap.add_argument("--no-conda", action="store_true")
args = ap.parse_args()

random.seed(11)
BASES = "ACGT"


def rand_seq(n):
    return "".join(random.choice(BASES) for _ in range(n))


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


def mutate(s, rate):
    return "".join(random.choice(BASES) if random.random() < rate else c for c in s)


if SB.exists():
    shutil.rmtree(SB)
SB.mkdir(parents=True)
for name in ("Snakefile", "workflow", "R", "envs"):
    (SB / name).symlink_to(ROOT / name)
(SB / "config").mkdir()
shutil.copy(ROOT / "config" / "config.yaml", SB / "config" / "config.yaml")

with open(ROOT / "config" / "samples.tsv") as fh:
    manifest = list(csv.DictReader(fh, delimiter="\t"))
manifest = manifest[:6]
with open(SB / "config" / "samples.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(manifest[0]), delimiter="\t")
    w.writeheader()
    w.writerows(manifest)
samples = [r["sample_id"] for r in manifest]

# ---------------------------------------------------------------- references
herv = rand_seq(6000)
viruses = {
    "NC_022518.1": ("Human endogenous retrovirus K113, complete genome", herv),
    "NC_900001.1": ("Fake phage A, complete genome", rand_seq(8000)),
    "NC_900002.1": ("Fake phage B, complete genome", rand_seq(8000)),
    "NC_900003.1": ("Fake virus C, complete genome", rand_seq(5000)),
}
host = rand_seq(30000) + herv + rand_seq(30000)

res = SB / "resources"
(res / "index").mkdir(parents=True)
(res / "hostile").mkdir()
with gzip.open(res / "refseq_viral.fna.gz", "wt") as fh:
    for acc, (desc, seq) in viruses.items():
        fh.write(f">{acc} {desc}\n{seq}\n")
(res / "refseq_viral.source.txt").write_text("url\tsynthetic\n")
with open(res / "refseq_viral_contigs.tsv", "w") as fh:
    fh.write("accession\tlength\tdescription\n")
    for acc, (desc, seq) in viruses.items():
        fh.write(f"{acc}\t{len(seq)}\t{desc}\n")
with open(res / "host.fa", "w") as fh:
    fh.write(f">host\n{host}\n")

subprocess.run(["bowtie2-build", "-q", str(res / "refseq_viral.fna.gz"),
                str(res / "index" / "refseq_viral")], check=True)
subprocess.run(["bowtie2-build", "-q", str(res / "host.fa"),
                str(res / "hostile" / "human-t2t-hla")], check=True)
(res / "hostile" / "human-t2t-hla.source.txt").write_text("index\tsynthetic\n")

# ---------------------------------------------------------------- reads
L, INSERT = 150, 350


def pair_from(seq, err=0.002):
    start = random.randint(0, len(seq) - INSERT)
    frag = seq[start:start + INSERT]
    if random.random() < 0.5:
        frag = revcomp(frag)
    return mutate(frag[:L], err), mutate(revcomp(frag)[:L], err)


planted = {}
fq = SB / "results" / "fastq_trimmed"
fq.mkdir(parents=True)
(SB / "results" / "qc" / "fastp").mkdir(parents=True)
for i, s in enumerate(samples):
    pairs = []
    n_host_herv = 40 + 5 * i
    for _ in range(300):
        pairs.append(pair_from(host[:30000]))
    for _ in range(n_host_herv):
        pairs.append(pair_from(herv))
    for acc in ("NC_900001.1", "NC_900002.1"):
        for _ in range(30 + 10 * (i % 3)):
            pairs.append(pair_from(viruses[acc][1]))
    if i % 2 == 0:
        for _ in range(25):
            # diverged reads: an aligner-sensitive case
            pairs.append(pair_from(viruses["NC_900003.1"][1], err=0.06))
    for _ in range(200):
        pairs.append((rand_seq(L), rand_seq(L)))
    random.shuffle(pairs)
    planted[s] = {"pairs": len(pairs), "herv_pairs": n_host_herv}
    with gzip.open(fq / f"{s}_1.fastq.gz", "wt") as f1, gzip.open(fq / f"{s}_2.fastq.gz", "wt") as f2:
        for j, (a, b) in enumerate(pairs):
            f1.write(f"@{s}.{j}\n{a}\n+\n{'I' * L}\n")
            f2.write(f"@{s}.{j}\n{b}\n+\n{'I' * L}\n")
    with open(SB / "results" / "qc" / "fastp" / f"{s}.json", "w") as fh:
        json.dump({"summary": {"before_filtering": {"total_reads": 2 * len(pairs)},
                               "after_filtering": {"total_reads": 2 * len(pairs)}}}, fh)

# ---------------------------------------------------------------- run
base = ["snakemake", "--cores", "4", "--rerun-triggers", "mtime", "--"]
# v1 tables first, exactly as the original rules make them, so the v2
# baseline reproduction check has something real to compare against.
subprocess.run(base + ["results/counts/virus_counts_long.tsv",
                       "results/qc/read_retention.tsv"], cwd=SB, check=True)
v2 = base[:-1]
if not args.no_conda:
    v2 += ["--use-conda", "--conda-prefix", str(ROOT / "build" / "conda-envs")]
v2 += ["--", "v2_all"]
subprocess.run(v2, cwd=SB, check=True)

# ---------------------------------------------------------------- checks
V2 = SB / "results" / "v2"
expected = [
    "host_depletion/host_depletion_per_sample.tsv",
    "host_depletion/reference_changes.tsv",
    "host_depletion/baseline_reproduction.tsv",
    "aligner_comparison/per_aligner.tsv",
    "aligner_comparison/per_sample.tsv",
    "aligner_comparison/breadth_per_reference.tsv",
    "aligner_comparison/agreement.tsv",
    "aligner_comparison/discordant_reads.tsv",
    "aligner_comparison/discordant_summary.tsv",
    "stats/dehost_bowtie2/permanova.tsv",
    "stats/baseline_bowtie2/permanova.tsv",
]
missing = [p for p in expected if not (V2 / p).exists()]
assert not missing, f"missing outputs: {missing}"

with open(V2 / "host_depletion" / "baseline_reproduction.tsv") as fh:
    repro = list(csv.DictReader(fh, delimiter="\t"))
assert repro and repro[0]["check"] == "PASS", repro

with open(V2 / "host_depletion" / "host_depletion_per_sample.tsv") as fh:
    rows = {r["sample_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
for s in samples:
    r = rows[s]
    assert r["input_matches_post_qc"] == "1", r
    assert int(r["host_reads_removed"]) > 0, r
    assert int(r["herv_k113_reads_before"]) > 0, r
    assert int(r["herv_k113_reads_after"]) == 0, r
    print(f"{s}: planted {planted[s]['herv_pairs']} HERV pairs; "
          f"HERV reads {r['herv_k113_reads_before']} -> {r['herv_k113_reads_after']}; "
          f"host removed {r['host_reads_removed']} reads ({r['host_pct_of_post_qc']}%)")

print((V2 / "aligner_comparison" / "per_aligner.tsv").read_text())
print((V2 / "aligner_comparison" / "discordant_summary.tsv").read_text())
print("smoke_v2: OK")
