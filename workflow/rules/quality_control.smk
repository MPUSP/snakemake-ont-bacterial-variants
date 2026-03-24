# --------------------------------------------------------------------------- #
# Raw read QC                                                                 #
# --------------------------------------------------------------------------- #
rule nanoplot_rawfastq:
    input:
        fastq=get_fastq,
    output:
        report=os.path.join(outdir, "qc/raw_reads/{sample}/report.html"),
        stats=os.path.join(outdir, "qc/raw_reads/{sample}/NanoStats.txt"),
        raw=os.path.join(outdir, "qc/raw_reads/{sample}/NanoPlot-data.tsv.gz"),
    log:
        os.path.join(outdir, "qc/raw_reads/logs/{sample}.log"),
    threads: config["nanoplot"]["threads"]
    params:
        extra="--raw",
    message:
        "--- Read QC, raw"
    wrapper:
        "https://raw.githubusercontent.com/MPUSP/mpusp-snakemake-wrappers/refs/heads/main/nanoplot"


# --------------------------------------------------------------------------- #
# Filtered read QC                                                            #
# --------------------------------------------------------------------------- #
rule nanoplot_filteredfastq:
    input:
        fastq=os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
    output:
        report=os.path.join(outdir, "qc/filtered_reads/{sample}/report.html"),
        stats=os.path.join(outdir, "qc/filtered_reads/{sample}/NanoStats.txt"),
        raw=os.path.join(outdir, "qc/filtered_reads/{sample}/NanoPlot-data.tsv.gz"),
    log:
        os.path.join(outdir, "qc/filtered_reads/logs/{sample}.log"),
    threads: config["nanoplot"]["threads"]
    params:
        extra="--raw",
    message:
        "--- Read QC, filtered"
    wrapper:
        "https://raw.githubusercontent.com/MPUSP/mpusp-snakemake-wrappers/refs/heads/main/nanoplot"


# --------------------------------------------------------------------------- #
# Aligned reads QC                                                            #
# --------------------------------------------------------------------------- #
rule nanoplot_aligned:
    input:
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
    output:
        report=os.path.join(outdir, "qc/aligned_reads/{sample}/report.html"),
        stats=os.path.join(outdir, "qc/aligned_reads/{sample}/NanoStats.txt"),
        raw=os.path.join(outdir, "qc/aligned_reads/{sample}/NanoPlot-data.tsv.gz"),
    log:
        os.path.join(outdir, "qc/aligned_reads/logs/{sample}.log"),
    threads: config["nanoplot"]["threads"]
    params:
        extra="--raw",
    message:
        "--- Read QC, aligned region"
    wrapper:
        "https://raw.githubusercontent.com/MPUSP/mpusp-snakemake-wrappers/refs/heads/main/nanoplot"


# --------------------------------------------------------------------------- #
# Merge QC reports using MultiQC                                              #
# --------------------------------------------------------------------------- #
rule multiqc:
    input:
        expand(rules.nanoplot_rawfastq.output.stats, sample=SAMPLES),
        expand(rules.nanoplot_filteredfastq.output.stats, sample=SAMPLES),
        expand(rules.nanoplot_aligned.output.stats, sample=SAMPLES),
    output:
        outdir=directory(os.path.join(outdir, "qc/multiqc")),
        report=os.path.join(outdir, "qc/multiqc/multiqc_report.html"),
        final=os.path.join(outdir, "variant_reports/MultiQC.html"),
    log:
        os.path.join(outdir, "qc/multiqc/multiqc.log"),
    conda:
        "../envs/multiqc.yml"
    threads: config["multiqc"]["threads"]
    params:
        qcdir=lambda x, output: os.path.dirname(str(output.outdir)),
    message:
        "--- MultiQC"
    shell:
        "multiqc "
        "--force "
        "--dirs "
        "--dirs-depth 2 "
        "--outdir {output.outdir} "
        "{params.qcdir} && "
        "cp {output.report} {output.final}"
