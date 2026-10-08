# results/v2

Host depletion and aligner comparison. The rules are in
`workflow/rules/v2.smk` and the design is described in the top-level README
under "Follow-up experiments".

## Status: rules written and tested, not yet run on the real data

There are no results in this directory yet. No number for either experiment
has been produced, so none is reported anywhere in this repository.

The rules were added from a cloud session whose network policy blocked the
three sources the run needs:

| Needed for | Host | Result |
|---|---|---|
| Hostile `human-t2t-hla` index | `objectstorage.uk-london-1.oraclecloud.com` | blocked (HTTP 403 at the proxy) |
| RefSeq viral release 237 | `ftp.ncbi.nlm.nih.gov` | blocked (HTTP 403 at the proxy) |
| SRA run lookup for `prefetch` | `trace.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov` | blocked (HTTP 403 at the proxy) |

The trimmed reads from the original run are not in the repository (they are
gitignored), so they could not be reused there either.

What was verified in that session:

- `pytest tests/` passes, including `tests/test_v2lib.py` for the summary and
  read-classification logic.
- `snakemake -n v2_all` builds the full DAG.
- `python tools/smoke_v2.py` ran every v2 rule through Snakemake with real
  Hostile 2.0.2, bowtie2 2.5.5, minimap2 2.31, samtools and R on planted
  synthetic reads, and the outputs matched what was planted: every planted
  HERV read was removed by depletion, and the re-run baseline reproduced the
  v1-style count table exactly. Those synthetic numbers live in `build/` and
  are not results.

## To produce the results

On a machine with the original `results/fastq_trimmed/` (or network access to
NCBI to regenerate them) and access to the Hostile index host:

```bash
conda env create -f envs/pipeline.yaml   # adds minimap2 to oral-virome
snakemake --cores 4 --use-conda v2_all
```

Resource needs for the Hostile step: Hostile's README states 4 GB of RAM for
short reads with bowtie2. The index size could not be read from Hostile's
manifest (that host was blocked); a bowtie2 index of a whole human assembly
is typically around 4 GB, and the compressed download sits alongside it while
unpacking, so leave roughly twice that free. Treat the disk figure as an
estimate.
