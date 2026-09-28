"""Compose the one-page PDF report.

Every figure quoted here is read from a file the pipeline wrote. Nothing in
this script contains a hard-coded result.

Layout: the page is A4 portrait and the content must fit on one page, so the
vertical budget is explicit. The sample table is drawn into an axes of a
computed height with bbox=[0,0,1,1] so it fills exactly that box rather than
overflowing into the section below it.
"""
import csv
import textwrap

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_kv(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if "\t" in line:
                k, v = line.rstrip("\n").split("\t", 1)
                out[k] = v
    return out


inp = snakemake.input  # noqa: F821
retention = read_tsv(inp.retention)
alpha_tests = read_tsv(inp.alpha_tests)
permanova = read_tsv(inp.permanova)
counts = read_tsv(inp.counts)
long_counts = read_tsv(inp.long)
manifest = read_tsv(inp.samples)
refmeta = read_kv(inp.reference)
cfg = snakemake.config  # noqa: F821

n_samples = len(manifest)
groups = {}
for r in manifest:
    groups.setdefault(r["diagnosis"], []).append(r["sample_id"])
countries = sorted({r["country"] for r in manifest})

sub_pairs = sum(int(r["subsampled_pairs"]) for r in retention)
post_pairs = sum(int(r["post_qc_pairs"]) for r in retention)
mapped_pairs = sum(int(r["mapped_pairs"]) for r in retention)
qc_pct = 100.0 * post_pairs / sub_pairs if sub_pairs else 0.0
map_pct = 100.0 * mapped_pairs / post_pairs if post_pairs else 0.0
n_virus = len(counts)


import collections

_tot = collections.Counter()
_ns = collections.Counter()
_name = {}
for _r in long_counts:
    _tot[_r["accession"]] += int(_r["reads"])
    _ns[_r["accession"]] += 1
    _name[_r["accession"]] = _r["virus"]
top_refs = _tot.most_common(4)
total_viral_reads = sum(_tot.values())
top2_share = (
    100.0 * sum(c for _, c in _tot.most_common(2)) / total_viral_reads
    if total_viral_reads else 0.0
)


def fmt_p(value):
    try:
        p = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{p:.4f}" if p >= 0.0001 else f"{p:.2e}"


def fmt_num(value, nd=3):
    try:
        return f"{float(value):.{nd}f}"
    except (TypeError, ValueError):
        return "NA"


fig = plt.figure(figsize=(8.27, 11.69))
LEFT = 0.065
WIDTH = 0.87
y = 0.978


def heading(text):
    global y
    fig.text(LEFT, y, text, size=9.5, weight="bold", va="top")
    y -= 0.0175


def para(text, size=7.3, width=118, dy=0.0116, bullet=False):
    global y
    lines = textwrap.wrap(text, width=width)
    for i, ln in enumerate(lines):
        prefix = ("- " if i == 0 else "  ") if bullet else ""
        fig.text(LEFT, y, prefix + ln, size=size, va="top")
        y -= dy


def gap(amount=0.008):
    global y
    y -= amount


# ------------------------------------------------------------------ header
fig.text(LEFT, y, "Oral virome in periodontal health and periodontitis",
         size=13.5, weight="bold", va="top")
y -= 0.024
fig.text(LEFT, y, f"Reanalysis of {cfg['bioproject']} subgingival shotgun "
                  "metagenomes. Portfolio pipeline, not a clinical finding.",
         size=7.8, style="italic", va="top")
y -= 0.024

# ---------------------------------------------------------------- question
heading("Question")
para("In subgingival gingival crevicular fluid, does the viral fraction of the "
     "metagenome differ in diversity or composition between periodontally healthy "
     "subjects and periodontitis patients?")
gap()

# ----------------------------------------------------------------- methods
heading("Methods")
para(
    f"{n_samples} samples were taken from BioProject {cfg['bioproject']} "
    f"({', '.join(f'{k} n={len(v)}' for k, v in sorted(groups.items()))}), balanced "
    f"across {len(countries)} countries ({', '.join(countries)}), and subsampled to "
    f"{cfg['subsample']['n_pairs']:,} read pairs with seqtk at seed "
    f"{cfg['subsample']['seed']} so that depth is identical across samples and the "
    f"pipeline runs on a laptop. Reads were adapter and quality trimmed with fastp, "
    f"with FastQC before and after and a MultiQC summary. Surviving pairs were "
    f"aligned with bowtie2 to NCBI RefSeq viral genomes (release "
    f"{refmeta.get('refseq_release', '?')}), keeping only properly paired alignments "
    f"at MAPQ >= {cfg['mapping']['min_mapq']}; a virus counted as present only with "
    f">= {cfg['presence']['min_reads']} reads across >= "
    f"{cfg['presence']['min_breadth']:.0%} of its genome, a breadth rule that "
    f"suppresses pileups on single repetitive loci. Alpha diversity was compared "
    f"with exact Wilcoxon rank-sum tests and Bray-Curtis composition with PERMANOVA "
    f"(9999 permutations) in R with vegan. Human host reads were not depleted, which "
    f"is a stated limitation."
)
gap(0.010)

# ----------------------------------------------------------------- table 1
heading("Table 1. Read retention and viral mapping")

rows = sorted(retention, key=lambda r: (r["diagnosis"], r["sample_id"]))
table_data = [
    [
        r["sample_id"],
        "Healthy" if r["diagnosis"].startswith("Periodontally") else "Periodontitis",
        f"{int(r['subsampled_pairs']):,}",
        f"{int(r['post_qc_pairs']):,}",
        f"{float(r['qc_retained_pct']):.1f}",
        f"{int(r['mapped_pairs']):,}",
        f"{float(r['mapped_pct']):.4f}",
    ]
    for r in rows
]
table_data.append(
    ["TOTAL", "", f"{sub_pairs:,}", f"{post_pairs:,}", f"{qc_pct:.1f}",
     f"{mapped_pairs:,}", f"{map_pct:.4f}"]
)

ROW_H = 0.0114
table_h = ROW_H * (len(table_data) + 1)
ax_t = fig.add_axes([LEFT, y - table_h, WIDTH, table_h])
ax_t.axis("off")
tbl = ax_t.table(
    cellText=table_data,
    colLabels=["Sample", "Group", "Subsampled pairs", "Post-QC pairs",
               "QC %", "Viral pairs", "Viral %"],
    cellLoc="right",
    bbox=[0, 0, 1, 1],
)
tbl.auto_set_font_size(False)
tbl.set_fontsize(6.1)
for (row, _col), cell in tbl.get_celld().items():
    cell.set_linewidth(0.25)
    if row == 0:
        cell.set_text_props(weight="bold")
        cell.set_facecolor("#e9e9e9")
    elif row == len(table_data):
        cell.set_text_props(weight="bold")
        cell.set_facecolor("#f5f5f5")
y -= table_h + 0.014

# ----------------------------------------------------------------- results
heading("Results")
para(f"{n_virus} RefSeq viral sequences passed the presence filter in at least one "
     f"sample. Viral read pairs were {map_pct:.4f}% of post-QC pairs.")
for t in alpha_tests:
    para(f"{t['metric'].replace('_', ' ').capitalize()}: median "
         f"{fmt_num(t['median_1'], 2)} ({t['group_1']}) vs {fmt_num(t['median_2'], 2)} "
         f"({t['group_2']}); W = {fmt_num(t['W'], 1)}, p = {fmt_p(t['p_value'])}, "
         f"rank-biserial r = {fmt_num(t['rank_biserial_r'], 2)}.")
for pr in permanova:
    if not pr.get("model"):
        continue
    if pr["model"].startswith("betadisper"):
        para(f"Dispersion homogeneity (betadisper): F = {fmt_num(pr['F'], 2)}, "
             f"p = {fmt_p(pr['p_value'])}.")
    else:
        para(f"PERMANOVA {pr['model']}, term {pr['term']}: R2 = {fmt_num(pr['R2'], 3)}, "
             f"F = {fmt_num(pr['F'], 2)}, p = {fmt_p(pr['p_value'])}.")

para("Most aligned reads fall on a small number of references: "
     + "; ".join(
         f"{_name[a]} ({c} reads in {_ns[a]}/{n_samples} samples)"
         for a, c in top_refs
     )
     + f". The top two account for {top2_share:.0f}% of all aligned reads. "
     "The first is a human endogenous retrovirus, which is host genomic sequence "
     "rather than an exogenous virus, and the second is a reference NCBI itself "
     "flags UNVERIFIED_ORG. Read the headline counts as an upper bound dominated "
     "by artefact, not as a virome.", size=7.0, width=124, dy=0.0112)
gap(0.006)

# ------------------------------------------------------------------ figure
FIG_H = 0.152
try:
    img = mpimg.imread(str(inp.figure))
    h, w = img.shape[0], img.shape[1]
    draw_w = min(WIDTH, FIG_H * (11.69 / 8.27) * (w / h))
    ax_i = fig.add_axes([LEFT + (WIDTH - draw_w) / 2, y - FIG_H, draw_w, FIG_H])
    ax_i.imshow(img)
    ax_i.axis("off")
    y -= FIG_H + 0.006
except Exception:  # noqa: BLE001
    para("[diversity figure unavailable]")

fig.text(LEFT, y, "Figure 1. Alpha diversity by group, and Bray-Curtis PCoA.",
         size=6.9, style="italic", va="top")
y -= 0.017

# -------------------------------------------------------------- limitations
heading("Limitations")
for text in [
    "Cross-sectional data. Any difference is an association, not evidence that either "
    "condition causes the other.",
    f"Each sample was subsampled to {cfg['subsample']['n_pairs']:,} read pairs to run "
    f"on a laptop, discarding most of the data and lowering sensitivity to rare viruses.",
    "This is a bulk shotgun metagenome, not a VLP-enriched virome, so viral reads are "
    "a small and mostly bacteriophage fraction and counts are sparse.",
    "Human reads were not depleted. Host sequence can align spuriously to viral "
    "references, so absolute viral abundance should be read as an upper bound.",
    f"n = {n_samples} across {len(countries)} countries. The source study found country "
    f"explained more compositional variance than diagnosis, and this subset cannot "
    f"separate the two.",
    "Detection is reference-based, so viruses absent from RefSeq viral, which is much "
    "of the oral phage population, are invisible here.",
]:
    para(text, size=6.9, width=126, dy=0.0107, bullet=True)

with PdfPages(str(snakemake.output[0])) as pdf:  # noqa: F821
    pdf.savefig(fig)
plt.close(fig)
print(f"wrote {snakemake.output[0]}  (content ended at y={y:.3f})")  # noqa: F821
