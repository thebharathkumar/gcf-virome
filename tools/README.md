# tools/

Helper scripts that are not part of the Snakemake workflow. They are kept in
the repository because they document how the dataset was chosen and verified,
which is part of the work.

## Dataset selection and provenance

| Script | Purpose |
|---|---|
| `ncbi.py` | Thin E-utilities client. Read-only. Every accession quoted anywhere in this repo was checked through it. |
| `survey_candidates.py` | Lists BioProjects matching oral and salivary virome queries. The first pass of the dataset search. |
| `verify_candidates.py` | For a shortlist of BioProjects, reports run counts, library strategy, data volume, per-sample attributes and any linked paper. This is what ruled out the alternatives. |
| `profile_target.py` | Full 80-sample profile of PRJNA1183294 and the PubMed search that matched it to its paper. |
| `fetch_abstracts.py` | Fetches abstracts for the candidate papers, used to confirm the paper matches the accession by study design. |
| `make_manifest.py` | Generates `config/samples.tsv` and `config/bioproject_all_samples.tsv` from live NCBI metadata. This is the only sanctioned way to create the manifest; it is never hand-edited. |

Why these are kept: the dataset was not picked by memory. `PRJNA692713`, for
example, looked ideal from its title ("microbiome and virome in oral cavity
squamous cell carcinoma, case control, shotgun") and turned out on inspection
to hold 9 runs totalling 31 spots. Re-running `verify_candidates.py` shows
that.

## Smoke tests

These exercise pipeline scripts outside Snakemake on synthetic input, so that
a bug in a late stage surfaces immediately rather than after a multi-hour run.
They write only to `build/`, never to `results/`, and their numbers are
meaningless.

| Script | Covers |
|---|---|
| `smoke_metadata.py` | `build_metadata.py` and `validate_metadata.py` against the real manifest |
| `smoke_diversity.py` | `R/diversity.R` on a synthetic count matrix |
| `smoke_report.py` | `build_count_matrix.py` and `make_report.py` |
| `smoke_v2.py` | Every `results/v2` rule end to end through Snakemake, with real Hostile, bowtie2, minimap2, samtools and R, on planted reads in a sandbox under `build/smoke_v2/` |

Run them with the pipeline environment active:

```bash
python tools/smoke_metadata.py
python tools/smoke_diversity.py
python tools/smoke_report.py
python tools/smoke_v2.py     # builds the Hostile env on first use
```

`smoke_metadata.py` is the only one whose output is meaningful, because it
runs on the real manifest. It exits non-zero if validation fails.
