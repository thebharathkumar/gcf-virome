# results/v2

Host depletion and aligner comparison. The rules are in
`workflow/rules/v2.smk` and the design is described in the top-level README
under "Follow-up experiments".

## Status: run on the real data 2026-10-08

Both experiments have been run to completion on all 16 samples. The headline
numbers and the limitations are in the top-level README under "Follow-up
experiments"; the files in this directory are the source for every one of
them.

| File | What it holds |
|---|---|
| `host_depletion/host_depletion_per_sample.tsv` | Host reads removed, viral pairs, HERV-K113 reads and references present, before and after, per sample plus a TOTAL row |
| `host_depletion/reference_changes.tsv` | Every (sample, reference) row whose reads or presence changed |
| `host_depletion/baseline_reproduction.tsv` | v1 count table against the v2 re-run of the original step |
| `aligner_comparison/per_aligner.tsv`, `per_sample.tsv` | Viral pairs, references present and wall time per aligner |
| `aligner_comparison/agreement.tsv` | Which aligner called each (sample, reference) |
| `aligner_comparison/discordant_reads.tsv`, `discordant_summary.tsv` | Per-read lookup of every discordant call in the other aligner's unfiltered BAM |
| `stats/baseline_bowtie2/`, `stats/dehost_bowtie2/` | `R/diversity.R` rerun on each count table |

The run took 13 hours 45 minutes wall clock on an Apple Silicon MacBook Air at
`--cores 4`, of which Hostile was 38,697.7 seconds (10.75 hours) per
`host_depletion/benchmark/*.tsv`.

Two notes on the environment, both of which affect what is in these files:

- Snakemake cannot read memory or IO counters on macOS, so `max_rss` and every
  related field in the `benchmark/` files is `NA`, and `max_rss_mb` in
  `aligner_comparison/per_sample.tsv` is `NA` for all 32 rows. The timing
  columns are real.
- `stats/dehost_bowtie2/permanova.tsv` is a single all-`NA` row. That is the
  honest output, not a failure: after depletion only 2 nonzero cells remain
  across 16 samples, so the Bray-Curtis dissimilarity is undefined for the
  other 14.

FASTQ and BAM files under this directory are gitignored; the summary tables,
counts, statistics and benchmarks are committed.

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
