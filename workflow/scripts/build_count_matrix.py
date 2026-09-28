"""Build the sample-by-virus count table from samtools coverage output.

Two tables come out of this:

  virus_count_matrix.tsv     raw aligned read counts per reference sequence
  virus_presence_matrix.tsv  0/1 after the pre-specified presence filter

The presence filter requires both a minimum read count and a minimum breadth
of coverage. Read count alone is not enough: in a bulk metagenome a single
repetitive or low-complexity locus can absorb a stack of reads and look like a
present virus. Requiring that the reads are spread over some fraction of the
genome removes most of that.
"""
import csv
import os
import re
from collections import defaultdict

coverage_files = snakemake.input.coverage  # noqa: F821
contigs_file = snakemake.input.contigs  # noqa: F821
samples_file = snakemake.input.samples  # noqa: F821
min_reads = int(snakemake.params.min_reads)  # noqa: F821
min_breadth = float(snakemake.params.min_breadth)  # noqa: F821

with open(samples_file) as fh:
    sample_ids = [r["sample_id"] for r in csv.DictReader(fh, delimiter="\t")]

descriptions = {}
with open(contigs_file) as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        descriptions[row["accession"]] = row["description"]


def virus_name(accession):
    """Trim the RefSeq description down to an organism-ish label."""
    desc = descriptions.get(accession, accession)
    # Drop trailing sequence-role wording so segments of one virus collapse
    # to the same label only when they genuinely share a name.
    desc = re.sub(
        r",?\s*(complete|partial)\s+(genome|sequence|cds)\.?$", "", desc, flags=re.I
    )
    return desc.strip() or accession


counts = defaultdict(dict)  # accession -> sample -> reads
breadth = defaultdict(dict)  # accession -> sample -> covered fraction

for path in coverage_files:
    sample = os.path.basename(path).replace(".coverage.tsv", "")
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            # samtools coverage writes a leading '#rname' header.
            acc = row.get("#rname") or row.get("rname")
            nreads = int(row["numreads"])
            if nreads == 0:
                continue
            counts[acc][sample] = nreads
            breadth[acc][sample] = float(row["coverage"]) / 100.0

observed = sorted(counts)
print(f"{len(observed)} reference sequences received at least one read")

presence = {}
for acc in observed:
    presence[acc] = {
        s: int(counts[acc].get(s, 0) >= min_reads and breadth[acc].get(s, 0.0) >= min_breadth)
        for s in sample_ids
    }

kept = [a for a in observed if any(presence[a].values())]
print(f"{len(kept)} pass the presence filter in at least one sample "
      f"(>= {min_reads} reads and >= {min_breadth:.1%} breadth)")


def write_matrix(path, rows, table):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["accession", "virus"] + sample_ids)
        for acc in rows:
            w.writerow([acc, virus_name(acc)] + [table[acc].get(s, 0) for s in sample_ids])


write_matrix(snakemake.output.counts, kept, counts)  # noqa: F821
write_matrix(snakemake.output.presence, kept, presence)  # noqa: F821

with open(snakemake.output.long, "w", newline="") as fh:  # noqa: F821
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample_id", "accession", "virus", "reads", "breadth", "present"])
    for acc in observed:
        for s in sample_ids:
            reads = counts[acc].get(s, 0)
            if reads == 0:
                continue
            w.writerow(
                [
                    s,
                    acc,
                    virus_name(acc),
                    reads,
                    f"{breadth[acc].get(s, 0.0):.6f}",
                    presence[acc][s],
                ]
            )

print(f"wrote {snakemake.output.counts}")  # noqa: F821
