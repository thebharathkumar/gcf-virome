"""Parse accession, length and description from the RefSeq viral FASTA."""
import gzip

fasta = snakemake.input[0]  # noqa: F821
out = snakemake.output[0]  # noqa: F821


def records(path):
    acc = desc = None
    length = 0
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc is not None:
                    yield acc, length, desc
                header = line[1:].strip()
                acc, _, desc = header.partition(" ")
                length = 0
            else:
                length += len(line.strip())
    if acc is not None:
        yield acc, length, desc


n = 0
with open(out, "w") as fh:
    fh.write("accession\tlength\tdescription\n")
    for acc, length, desc in records(fasta):
        fh.write(f"{acc}\t{length}\t{desc}\n")
        n += 1

print(f"parsed {n} reference sequences from {fasta}")
