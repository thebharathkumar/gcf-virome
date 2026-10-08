"""results/v2: host depletion and aligner comparison.

Two experiments layered on the original run without touching it. Every output
here lives under results/v2/ (or resources/ for reference indexes), and every
filter is read from the same config keys the original rules use, so the
comparison is like for like.

  1. Host depletion. Hostile (bowtie2 mode, human-t2t-hla index) removes human
     read pairs from the trimmed reads, then the original bowtie2 alignment is
     repeated on what is left.
  2. Aligner comparison. The host-depleted reads are aligned with bowtie2
     (original settings) and minimap2 -ax sr, then filtered and counted
     identically.

A third alignment, "baseline_bowtie2", repeats the original bowtie2 step on
the non-depleted reads inside the v2 tree. It is the "before" side of
experiment 1, and it is checked against the committed v1 count table so any
drift between machines or tool builds is reported rather than hidden.

Run with:  snakemake --cores 4 --use-conda v2_all

--use-conda matters only for the two Hostile rules, which need their own
environment (envs/hostile.yaml). Every other rule runs in the main
oral-virome environment exactly as the original pipeline does.
"""

V2 = config["v2"]["outdir"]
HOSTILE_INDEX = config["v2"]["hostile"]["index_name"]
HOSTILE_DIR = "resources/hostile"
HOSTILE_PREFIX = f"{HOSTILE_DIR}/{HOSTILE_INDEX}"
BT2_EXT = [".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"]

# run name -> (input reads, aligner)
V2_RUNS = {
    "baseline_bowtie2": ("trimmed", "bowtie2"),
    "dehost_bowtie2": ("dehost", "bowtie2"),
    "dehost_minimap2": ("dehost", "minimap2"),
}
# Diversity statistics are rerun on the depleted bowtie2 counts, and on the
# v2 baseline so before and after come from the same software build.
V2_DIVERSITY_RUNS = ["baseline_bowtie2", "dehost_bowtie2"]


wildcard_constraints:
    run="|".join(V2_RUNS),


def v2_reads(wc):
    source = V2_RUNS[wc.run][0]
    if source == "trimmed":
        d = "results/fastq_trimmed"
    else:
        d = f"{V2}/host_depletion/fastq"
    return {"r1": f"{d}/{wc.sample}_1.fastq.gz", "r2": f"{d}/{wc.sample}_2.fastq.gz"}


def v2_index(wc):
    if V2_RUNS[wc.run][1] == "bowtie2":
        return multiext("resources/index/refseq_viral", *BT2_EXT)
    return "resources/index/refseq_viral.sr.mmi"


rule v2_all:
    input:
        f"{V2}/host_depletion/host_depletion_per_sample.tsv",
        f"{V2}/host_depletion/baseline_reproduction.tsv",
        f"{V2}/aligner_comparison/per_aligner.tsv",
        f"{V2}/aligner_comparison/discordant_reads.tsv",
        expand(f"{V2}/stats/{{run}}/permanova.tsv", run=V2_DIVERSITY_RUNS),
        expand(f"{V2}/stats/{{run}}/alpha_group_tests.tsv", run=V2_DIVERSITY_RUNS),


# ---------------------------------------------------------------------------
# Experiment 1: host depletion
# ---------------------------------------------------------------------------

rule hostile_index:
    """Fetch the prebuilt Hostile index. Hostile verifies the SHA-256 itself."""
    output:
        multiext(HOSTILE_PREFIX, *BT2_EXT),
        meta=f"{HOSTILE_DIR}/{HOSTILE_INDEX}.source.txt",
    params:
        name=HOSTILE_INDEX,
        cache=HOSTILE_DIR,
    conda:
        "../../envs/hostile.yaml"
    log:
        "logs/v2/hostile_index.log",
    shell:
        """
        mkdir -p {params.cache}
        HOSTILE_CACHE_DIR={params.cache} hostile index fetch \
            --name {params.name} --bowtie2 > {log} 2>&1
        printf 'index\\t%s\\nhostile_version\\t%s\\nretrieved\\t%s\\n' \
            "{params.name}" "$(hostile --version 2>&1 | tail -1)" \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > {output.meta}
        """


rule hostile_clean:
    """Drop every pair in which either mate aligns to the human index.

    Hostile's paired bowtie2 mode keeps a pair only when both mates are
    unaligned (samtools view -f 12). Reads are not renamed, so read names stay
    traceable to the original FASTQ.
    """
    input:
        r1="results/fastq_trimmed/{sample}_1.fastq.gz",
        r2="results/fastq_trimmed/{sample}_2.fastq.gz",
        index=multiext(HOSTILE_PREFIX, *BT2_EXT),
    output:
        r1=f"{V2}/host_depletion/fastq/{{sample}}_1.fastq.gz",
        r2=f"{V2}/host_depletion/fastq/{{sample}}_2.fastq.gz",
        json=f"{V2}/host_depletion/hostile/{{sample}}.json",
    params:
        prefix=HOSTILE_PREFIX,
        tmp=f"{V2}/host_depletion/tmp/{{sample}}",
    threads: 4
    conda:
        "../../envs/hostile.yaml"
    benchmark:
        f"{V2}/host_depletion/benchmark/{{sample}}.tsv"
    log:
        "logs/v2/hostile_clean/{sample}.log",
    shell:
        """
        rm -rf {params.tmp} && mkdir -p {params.tmp}
        hostile clean \
            --fastq1 {input.r1} --fastq2 {input.r2} \
            --aligner bowtie2 --index {params.prefix} \
            --threads {threads} --airplane --force \
            --output {params.tmp} > {output.json} 2> {log}
        mv {params.tmp}/{wildcards.sample}_1.clean_1.fastq.gz {output.r1}
        mv {params.tmp}/{wildcards.sample}_2.clean_2.fastq.gz {output.r2}
        rm -rf {params.tmp}
        """


# ---------------------------------------------------------------------------
# Alignment, shared by both experiments
# ---------------------------------------------------------------------------

rule minimap2_index:
    input:
        "resources/refseq_viral.fna.gz",
    output:
        "resources/index/refseq_viral.sr.mmi",
    threads: 4
    benchmark:
        f"{V2}/aligner_comparison/benchmark/minimap2_index.tsv"
    log:
        "logs/v2/minimap2_index.log",
    shell:
        "mkdir -p resources/index && "
        "minimap2 -x sr -t {threads} -d {output} {input} > {log} 2>&1"


rule v2_align:
    """Align and keep every aligned record, unfiltered.

    The unfiltered BAM is kept so a read counted by one aligner can be looked
    up in the other aligner's output. Both aligners drop pairs where neither
    mate aligns (bowtie2 --no-unal, minimap2 --sam-hit-only). The benchmark
    covers alignment plus the identical samtools sort, with the same thread
    count for both aligners.
    """
    input:
        unpack(v2_reads),
        index=v2_index,
    output:
        bam=f"{V2}/align/{{run}}/raw/{{sample}}.bam",
        bai=f"{V2}/align/{{run}}/raw/{{sample}}.bam.bai",
    params:
        aligner=lambda wc: V2_RUNS[wc.run][1],
        bt2_prefix="resources/index/refseq_viral",
        preset=config["v2"]["minimap2_preset"],
    threads: 4
    benchmark:
        f"{V2}/align/{{run}}/benchmark/{{sample}}.tsv"
    log:
        "logs/v2/align/{run}/{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.bam})
        if [ "{params.aligner}" = "bowtie2" ]; then
            bowtie2 -x {params.bt2_prefix} -1 {input.r1} -2 {input.r2} \
                --no-unal --threads {threads} 2> {log} \
              | samtools sort -@ {threads} -o {output.bam} -
        else
            minimap2 -ax {params.preset} --sam-hit-only -t {threads} \
                {input.index} {input.r1} {input.r2} 2> {log} \
              | samtools sort -@ {threads} -o {output.bam} -
        fi
        samtools index {output.bam}
        """


rule v2_filter:
    """The original filter: MAPQ >= min_mapq and properly paired."""
    input:
        f"{V2}/align/{{run}}/raw/{{sample}}.bam",
    output:
        bam=f"{V2}/align/{{run}}/bam/{{sample}}.bam",
        bai=f"{V2}/align/{{run}}/bam/{{sample}}.bam.bai",
        stats=f"{V2}/align/{{run}}/bam/{{sample}}.flagstat",
    params:
        mapq=MIN_MAPQ,
    log:
        "logs/v2/filter/{run}/{sample}.log",
    shell:
        """
        mkdir -p $(dirname {output.bam})
        samtools view -b -q {params.mapq} -f 2 -o {output.bam} {input} 2> {log}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stats}
        """


rule v2_coverage:
    input:
        bam=f"{V2}/align/{{run}}/bam/{{sample}}.bam",
        bai=f"{V2}/align/{{run}}/bam/{{sample}}.bam.bai",
    output:
        f"{V2}/align/{{run}}/coverage/{{sample}}.coverage.tsv",
    log:
        "logs/v2/coverage/{run}/{sample}.log",
    shell:
        "samtools coverage {input.bam} > {output} 2> {log}"


rule v2_count_matrix:
    """Same script and same thresholds as the original count_matrix rule."""
    input:
        coverage=expand(f"{V2}/align/{{{{run}}}}/coverage/{{s}}.coverage.tsv", s=SAMPLES),
        contigs="resources/refseq_viral_contigs.tsv",
        samples=config["samples_tsv"],
    output:
        counts=f"{V2}/counts/{{run}}/virus_count_matrix.tsv",
        presence=f"{V2}/counts/{{run}}/virus_presence_matrix.tsv",
        long=f"{V2}/counts/{{run}}/virus_counts_long.tsv",
    params:
        min_reads=config["presence"]["min_reads"],
        min_breadth=config["presence"]["min_breadth"],
    log:
        "logs/v2/count_matrix/{run}.log",
    script:
        "../scripts/build_count_matrix.py"


rule v2_diversity:
    """The original R script, unchanged, pointed at v2 count tables."""
    input:
        counts=f"{V2}/counts/{{run}}/virus_count_matrix.tsv",
        presence=f"{V2}/counts/{{run}}/virus_presence_matrix.tsv",
        samples=config["samples_tsv"],
    output:
        alpha=f"{V2}/stats/{{run}}/alpha_diversity.tsv",
        alpha_tests=f"{V2}/stats/{{run}}/alpha_group_tests.tsv",
        permanova=f"{V2}/stats/{{run}}/permanova.tsv",
        pcoa=f"{V2}/stats/{{run}}/pcoa_coordinates.tsv",
        figure=f"{V2}/stats/{{run}}/diversity.png",
        session=f"{V2}/stats/{{run}}/r_session_info.txt",
    params:
        group=config["group_column"],
    log:
        "logs/v2/diversity/{run}.log",
    shell:
        "mkdir -p $(dirname {output.alpha}) && "
        "Rscript R/diversity.R {input.counts} {input.presence} {input.samples} "
        "{params.group} "
        "{output.alpha} {output.alpha_tests} {output.permanova} {output.pcoa} "
        "{output.figure} {output.session} > {log} 2>&1"


# ---------------------------------------------------------------------------
# Experiment 1 summary
# ---------------------------------------------------------------------------

rule v2_baseline_reproduction:
    """Does the v2 baseline reproduce the committed v1 count table exactly?"""
    input:
        v1_long="results/counts/virus_counts_long.tsv",
        v1_retention="results/qc/read_retention.tsv",
        v2_long=f"{V2}/counts/baseline_bowtie2/virus_counts_long.tsv",
        v2_flagstat=expand(f"{V2}/align/baseline_bowtie2/bam/{{s}}.flagstat", s=SAMPLES),
    output:
        f"{V2}/host_depletion/baseline_reproduction.tsv",
    log:
        "logs/v2/baseline_reproduction.log",
    script:
        "../scripts/v2_baseline_reproduction.py"


rule v2_host_depletion_summary:
    input:
        fastp=expand("results/qc/fastp/{s}.json", s=SAMPLES),
        hostile=expand(f"{V2}/host_depletion/hostile/{{s}}.json", s=SAMPLES),
        flagstat_before=expand(f"{V2}/align/baseline_bowtie2/bam/{{s}}.flagstat", s=SAMPLES),
        flagstat_after=expand(f"{V2}/align/dehost_bowtie2/bam/{{s}}.flagstat", s=SAMPLES),
        long_before=f"{V2}/counts/baseline_bowtie2/virus_counts_long.tsv",
        long_after=f"{V2}/counts/dehost_bowtie2/virus_counts_long.tsv",
        samples=config["samples_tsv"],
    output:
        per_sample=f"{V2}/host_depletion/host_depletion_per_sample.tsv",
        reference_changes=f"{V2}/host_depletion/reference_changes.tsv",
    params:
        herv=config["v2"]["herv_k113_accession"],
    log:
        "logs/v2/host_depletion_summary.log",
    script:
        "../scripts/v2_host_depletion_summary.py"


# ---------------------------------------------------------------------------
# Experiment 2 summary
# ---------------------------------------------------------------------------

rule v2_aligner_comparison:
    input:
        flagstat_bt2=expand(f"{V2}/align/dehost_bowtie2/bam/{{s}}.flagstat", s=SAMPLES),
        flagstat_mm2=expand(f"{V2}/align/dehost_minimap2/bam/{{s}}.flagstat", s=SAMPLES),
        bench_bt2=expand(f"{V2}/align/dehost_bowtie2/benchmark/{{s}}.tsv", s=SAMPLES),
        bench_mm2=expand(f"{V2}/align/dehost_minimap2/benchmark/{{s}}.tsv", s=SAMPLES),
        long_bt2=f"{V2}/counts/dehost_bowtie2/virus_counts_long.tsv",
        long_mm2=f"{V2}/counts/dehost_minimap2/virus_counts_long.tsv",
        contigs="resources/refseq_viral_contigs.tsv",
        samples=config["samples_tsv"],
    output:
        per_aligner=f"{V2}/aligner_comparison/per_aligner.tsv",
        per_sample=f"{V2}/aligner_comparison/per_sample.tsv",
        breadth=f"{V2}/aligner_comparison/breadth_per_reference.tsv",
        agreement=f"{V2}/aligner_comparison/agreement.tsv",
    log:
        "logs/v2/aligner_comparison.log",
    script:
        "../scripts/v2_aligner_comparison.py"


rule v2_explain_discordance:
    """Look up, read by read, why the aligners disagree.

    For every (sample, reference) called present by one aligner and not the
    other, each read the detecting aligner counted is looked up in the other
    aligner's unfiltered BAM and classified (unaligned, other reference, low
    MAPQ, not properly paired, or passing).
    """
    input:
        agreement=f"{V2}/aligner_comparison/agreement.tsv",
        raw_bt2=expand(f"{V2}/align/dehost_bowtie2/raw/{{s}}.bam", s=SAMPLES),
        raw_mm2=expand(f"{V2}/align/dehost_minimap2/raw/{{s}}.bam", s=SAMPLES),
        bam_bt2=expand(f"{V2}/align/dehost_bowtie2/bam/{{s}}.bam", s=SAMPLES),
        bam_mm2=expand(f"{V2}/align/dehost_minimap2/bam/{{s}}.bam", s=SAMPLES),
        contigs="resources/refseq_viral_contigs.tsv",
    output:
        reads=f"{V2}/aligner_comparison/discordant_reads.tsv",
        summary=f"{V2}/aligner_comparison/discordant_summary.tsv",
    params:
        dir=f"{V2}/align",
        mapq=MIN_MAPQ,
    log:
        "logs/v2/explain_discordance.log",
    script:
        "../scripts/v2_explain_discordance.py"
