"""Oral virome pipeline for PRJNA1183294.

Run with:  snakemake --cores 4

Disk note: the intermediate uncompressed FASTQ from fasterq-dump is large
(roughly 6 GB per sample). Those files are marked temp() and Snakemake deletes
them as soon as the subsampled copy exists, so peak usage stays near 10 GB
rather than 100 GB.
"""
import csv

configfile: "config/config.yaml"


def read_samples(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


SAMPLE_ROWS = read_samples(config["samples_tsv"])
SAMPLES = [r["sample_id"] for r in SAMPLE_ROWS]
RUN = {r["sample_id"]: r["run_accession"] for r in SAMPLE_ROWS}

N_PAIRS = config["subsample"]["n_pairs"]
SEED = config["subsample"]["seed"]
MIN_MAPQ = config["mapping"]["min_mapq"]


rule all:
    input:
        "results/qc/multiqc_report.html",
        "results/qc/read_retention.tsv",
        "results/counts/virus_count_matrix.tsv",
        "results/counts/virus_presence_matrix.tsv",
        "results/metadata/validation_report.txt",
        "results/stats/alpha_diversity.tsv",
        "results/stats/alpha_group_tests.tsv",
        "results/stats/permanova.tsv",
        "results/stats/pcoa_coordinates.tsv",
        "results/figures/diversity.png",
        "report/report.pdf",
        "results/summary.md",


# ---------------------------------------------------------------------------
# 1. Retrieval
# ---------------------------------------------------------------------------

rule prefetch:
    output:
        temp(directory("results/sra/{sample}")),
    params:
        run=lambda w: RUN[w.sample],
    log:
        "logs/prefetch/{sample}.log",
    shell:
        "prefetch {params.run} --max-size 50G -O {output} > {log} 2>&1"


rule fasterq_dump:
    input:
        "results/sra/{sample}",
    output:
        r1=temp("results/fastq_raw/{sample}_1.fastq"),
        r2=temp("results/fastq_raw/{sample}_2.fastq"),
    params:
        run=lambda w: RUN[w.sample],
    threads: 2
    log:
        "logs/fasterq_dump/{sample}.log",
    shell:
        """
        mkdir -p results/fastq_raw results/tmp
        fasterq-dump {input}/{params.run} \
            --split-files --threads {threads} \
            --temp results/tmp -O results/fastq_raw > {log} 2>&1
        mv results/fastq_raw/{params.run}_1.fastq {output.r1}
        mv results/fastq_raw/{params.run}_2.fastq {output.r2}
        """


# ---------------------------------------------------------------------------
# 2. Subsampling to a fixed depth with a fixed seed
# ---------------------------------------------------------------------------

rule subsample:
    """seqtk with one seed for both mates keeps the pairs in register."""
    input:
        r1="results/fastq_raw/{sample}_1.fastq",
        r2="results/fastq_raw/{sample}_2.fastq",
    output:
        r1="results/fastq_sub/{sample}_1.fastq.gz",
        r2="results/fastq_sub/{sample}_2.fastq.gz",
    params:
        n=N_PAIRS,
        seed=SEED,
    threads: 2
    log:
        "logs/subsample/{sample}.log",
    shell:
        """
        seqtk sample -s{params.seed} {input.r1} {params.n} | pigz -p {threads} > {output.r1} 2> {log}
        seqtk sample -s{params.seed} {input.r2} {params.n} | pigz -p {threads} > {output.r2} 2>> {log}
        """


# ---------------------------------------------------------------------------
# 3. QC before and after trimming
# ---------------------------------------------------------------------------

rule fastqc_raw:
    input:
        r1="results/fastq_sub/{sample}_1.fastq.gz",
        r2="results/fastq_sub/{sample}_2.fastq.gz",
    output:
        "results/qc/fastqc_raw/{sample}_1_fastqc.zip",
        "results/qc/fastqc_raw/{sample}_2_fastqc.zip",
    threads: 2
    log:
        "logs/fastqc_raw/{sample}.log",
    shell:
        "mkdir -p results/qc/fastqc_raw && "
        "fastqc -t {threads} -o results/qc/fastqc_raw {input.r1} {input.r2} > {log} 2>&1"


rule fastp:
    input:
        r1="results/fastq_sub/{sample}_1.fastq.gz",
        r2="results/fastq_sub/{sample}_2.fastq.gz",
    output:
        r1="results/fastq_trimmed/{sample}_1.fastq.gz",
        r2="results/fastq_trimmed/{sample}_2.fastq.gz",
        json="results/qc/fastp/{sample}.json",
        html="results/qc/fastp/{sample}.html",
    threads: 4
    log:
        "logs/fastp/{sample}.log",
    shell:
        """
        mkdir -p results/fastq_trimmed results/qc/fastp
        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --detect_adapter_for_pe \
            --cut_front --cut_tail --cut_mean_quality 20 \
            --qualified_quality_phred 20 --unqualified_percent_limit 30 \
            --length_required 50 \
            --thread {threads} \
            --json {output.json} --html {output.html} > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        r1="results/fastq_trimmed/{sample}_1.fastq.gz",
        r2="results/fastq_trimmed/{sample}_2.fastq.gz",
    output:
        "results/qc/fastqc_trimmed/{sample}_1_fastqc.zip",
        "results/qc/fastqc_trimmed/{sample}_2_fastqc.zip",
    threads: 2
    log:
        "logs/fastqc_trimmed/{sample}.log",
    shell:
        "mkdir -p results/qc/fastqc_trimmed && "
        "fastqc -t {threads} -o results/qc/fastqc_trimmed {input.r1} {input.r2} > {log} 2>&1"


rule multiqc:
    input:
        expand("results/qc/fastqc_raw/{s}_{r}_fastqc.zip", s=SAMPLES, r=[1, 2]),
        expand("results/qc/fastqc_trimmed/{s}_{r}_fastqc.zip", s=SAMPLES, r=[1, 2]),
        expand("results/qc/fastp/{s}.json", s=SAMPLES),
    output:
        "results/qc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    shell:
        "multiqc results/qc -o results/qc -n multiqc_report.html --force > {log} 2>&1"


# ---------------------------------------------------------------------------
# 4. Viral reference and alignment
# ---------------------------------------------------------------------------

rule download_viral_reference:
    output:
        fna="resources/refseq_viral.fna.gz",
        meta="resources/refseq_viral.source.txt",
    params:
        url=config["reference"]["viral_fna_url"],
        release=config["reference"]["refseq_release"],
    log:
        "logs/download_reference.log",
    shell:
        """
        mkdir -p resources
        curl -sSL --fail -o {output.fna} {params.url} 2> {log}
        printf 'url\\t%s\\nrefseq_release\\t%s\\nsha256\\t%s\\nbytes\\t%s\\nretrieved\\t%s\\n' \
            "{params.url}" "{params.release}" \
            "$(shasum -a 256 {output.fna} | cut -d' ' -f1)" \
            "$(wc -c < {output.fna} | tr -d ' ')" \
            "$(date -u +%Y-%m-%dT%H:%M:%SZ)" > {output.meta}
        """


rule bowtie2_index:
    input:
        "resources/refseq_viral.fna.gz",
    output:
        multiext(
            "resources/index/refseq_viral",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    threads: 4
    log:
        "logs/bowtie2_index.log",
    shell:
        "mkdir -p resources/index && "
        "bowtie2-build --threads {threads} {input} resources/index/refseq_viral > {log} 2>&1"


rule reference_lengths:
    """Accession -> (length, description), parsed from the FASTA headers."""
    input:
        "resources/refseq_viral.fna.gz",
    output:
        "resources/refseq_viral_contigs.tsv",
    log:
        "logs/reference_lengths.log",
    script:
        "workflow/scripts/reference_lengths.py"


rule bowtie2_map:
    input:
        r1="results/fastq_trimmed/{sample}_1.fastq.gz",
        r2="results/fastq_trimmed/{sample}_2.fastq.gz",
        index=multiext(
            "resources/index/refseq_viral",
            ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
        ),
    output:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai",
        stats="results/bam/{sample}.flagstat",
    params:
        prefix="resources/index/refseq_viral",
        mapq=MIN_MAPQ,
    threads: 4
    log:
        "logs/bowtie2_map/{sample}.log",
    shell:
        """
        mkdir -p results/bam
        bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} \
            --no-unal --threads {threads} 2> {log} \
          | samtools view -b -q {params.mapq} -f 2 - \
          | samtools sort -@ {threads} -o {output.bam} -
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stats}
        """


rule coverage:
    input:
        bam="results/bam/{sample}.bam",
        bai="results/bam/{sample}.bam.bai",
    output:
        "results/coverage/{sample}.coverage.tsv",
    log:
        "logs/coverage/{sample}.log",
    shell:
        "mkdir -p results/coverage && "
        "samtools coverage {input.bam} > {output} 2> {log}"


# ---------------------------------------------------------------------------
# 5. Count table and read retention
# ---------------------------------------------------------------------------

rule count_matrix:
    input:
        coverage=expand("results/coverage/{s}.coverage.tsv", s=SAMPLES),
        contigs="resources/refseq_viral_contigs.tsv",
        samples=config["samples_tsv"],
    output:
        counts="results/counts/virus_count_matrix.tsv",
        presence="results/counts/virus_presence_matrix.tsv",
        long="results/counts/virus_counts_long.tsv",
    params:
        min_reads=config["presence"]["min_reads"],
        min_breadth=config["presence"]["min_breadth"],
    log:
        "logs/count_matrix.log",
    script:
        "workflow/scripts/build_count_matrix.py"


rule read_retention:
    input:
        fastp=expand("results/qc/fastp/{s}.json", s=SAMPLES),
        flagstat=expand("results/bam/{s}.flagstat", s=SAMPLES),
        samples=config["samples_tsv"],
    output:
        "results/qc/read_retention.tsv",
    log:
        "logs/read_retention.log",
    script:
        "workflow/scripts/read_retention.py"


# ---------------------------------------------------------------------------
# 6. Metadata deliverables and validation
# ---------------------------------------------------------------------------

rule metadata_exports:
    input:
        samples=config["samples_tsv"],
    output:
        dictionary="results/metadata/redcap_data_dictionary.csv",
        records="results/metadata/redcap_records.csv",
        biosample="results/metadata/biosample_MIMS_human_oral.tsv",
        sra="results/metadata/sra_metadata.tsv",
    log:
        "logs/metadata_exports.log",
    script:
        "workflow/scripts/build_metadata.py"


rule validate_metadata:
    input:
        dictionary="results/metadata/redcap_data_dictionary.csv",
        records="results/metadata/redcap_records.csv",
        biosample="results/metadata/biosample_MIMS_human_oral.tsv",
        sra="results/metadata/sra_metadata.tsv",
        samples=config["samples_tsv"],
    output:
        "results/metadata/validation_report.txt",
    log:
        "logs/validate_metadata.log",
    shell:
        "mkdir -p results/metadata && "
        "python workflow/scripts/validate_metadata.py "
        "--dictionary {input.dictionary} --records {input.records} "
        "--biosample {input.biosample} --sra {input.sra} "
        "--manifest {input.samples} "
        "--report {output} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 7. Statistics in R
# ---------------------------------------------------------------------------

rule diversity:
    input:
        counts="results/counts/virus_count_matrix.tsv",
        presence="results/counts/virus_presence_matrix.tsv",
        samples=config["samples_tsv"],
    output:
        alpha="results/stats/alpha_diversity.tsv",
        alpha_tests="results/stats/alpha_group_tests.tsv",
        permanova="results/stats/permanova.tsv",
        pcoa="results/stats/pcoa_coordinates.tsv",
        figure="results/figures/diversity.png",
        session="results/stats/r_session_info.txt",
    params:
        group=config["group_column"],
    log:
        "logs/diversity.log",
    shell:
        "mkdir -p results/stats results/figures && "
        "Rscript R/diversity.R {input.counts} {input.presence} {input.samples} "
        "{params.group} "
        "{output.alpha} {output.alpha_tests} {output.permanova} {output.pcoa} "
        "{output.figure} {output.session} > {log} 2>&1"


# ---------------------------------------------------------------------------
# 8. One-page report
# ---------------------------------------------------------------------------

rule report:
    input:
        retention="results/qc/read_retention.tsv",
        alpha="results/stats/alpha_diversity.tsv",
        alpha_tests="results/stats/alpha_group_tests.tsv",
        permanova="results/stats/permanova.tsv",
        pcoa="results/stats/pcoa_coordinates.tsv",
        counts="results/counts/virus_count_matrix.tsv",
        presence="results/counts/virus_presence_matrix.tsv",
        long="results/counts/virus_counts_long.tsv",
        figure="results/figures/diversity.png",
        reference="resources/refseq_viral.source.txt",
        samples=config["samples_tsv"],
    output:
        "report/report.pdf",
    log:
        "logs/report.log",
    script:
        "workflow/scripts/make_report.py"


rule summary:
    """Regenerate the results table in README.md from the pipeline outputs."""
    input:
        retention="results/qc/read_retention.tsv",
        alpha="results/stats/alpha_diversity.tsv",
        alpha_tests="results/stats/alpha_group_tests.tsv",
        permanova="results/stats/permanova.tsv",
        counts="results/counts/virus_count_matrix.tsv",
        long="results/counts/virus_counts_long.tsv",
        samples=config["samples_tsv"],
    output:
        summary="results/summary.md",
    params:
        # README.md is patched as a side effect rather than declared as an
        # output. Snakemake deletes outputs when a job fails, and losing the
        # README to a transient error is not an acceptable failure mode.
        readme="README.md",
    log:
        "logs/summary.log",
    script:
        "workflow/scripts/write_summary.py"
