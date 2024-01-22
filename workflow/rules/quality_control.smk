# --------------------------------------------------------------------------- #
# Raw read QC                                                                 #
# --------------------------------------------------------------------------- #
rule nanoplot_rawfastq:
    """QC on raw nanopore reads"""
    message:
        "--- Read QC, raw"
    input:
        get_fastq
    output:
        outdir = directory(os.path.join(outdir, "qc/raw_reads/{sample}")),
        stats = os.path.join(outdir, "qc/raw_reads/{sample}/NanoStats.txt"),
        raw = os.path.join(outdir, "qc/raw_reads/{sample}/NanoPlot-data.tsv.gz")
    log:
        os.path.join(outdir, "qc/raw_reads/logs/{sample}.log")
    threads:
        config["nanoplot"]["threads"]
    conda:
        "../envs/nanoplot.yml"
    shell:
        "NanoPlot "
            "--threads {threads} "
            "--raw "
            "--fastq {input} "
            "--outdir {output.outdir} "
            "> {log}"

# --------------------------------------------------------------------------- #
# Filtered read QC                                                            #
# --------------------------------------------------------------------------- #
rule nanoplot_filteredfastq:
    """QC on filtered nanopore reads"""
    message:
        "--- Read QC, filtered"
    input:
        os.path.join(outdir, "filtered_reads/{sample}.fastq.gz")
    output:
        outdir = directory(os.path.join(outdir, "qc/filtered_reads/{sample}")),
        stats = os.path.join(outdir, "qc/filtered_reads/{sample}/NanoStats.txt"),
        raw = os.path.join(outdir, "qc/filtered_reads/{sample}/NanoPlot-data.tsv.gz")
    log:
        os.path.join(outdir, "qc/filtered_reads/logs/{sample}.log")
    threads:
        config["nanoplot"]["threads"]
    conda:
        "../envs/nanoplot.yml"
    shell:
        "NanoPlot "
            "--threads {threads} "
            "--raw "
            "--fastq {input} "
            "--outdir {output.outdir} "
            "> {log}"

# --------------------------------------------------------------------------- #
# Aligned reads QC                                                            #
# --------------------------------------------------------------------------- #
rule nanoplot_aligned:
    """QC on aligned read regions"""
    message:
        "--- Read QC, aligned region"
    input:
        os.path.join(outdir, "mapping/{sample}.bam")
    output:
        outdir = directory(os.path.join(outdir, "qc/aligned_reads/{sample}")),
        stats = os.path.join(outdir, "qc/aligned_reads/{sample}/NanoStats.txt"),
        raw = os.path.join(outdir, "qc/aligned_reads/{sample}/NanoPlot-data.tsv.gz")
    log:
        os.path.join(outdir, "qc/aligned_reads/logs/{sample}.log")
    threads:
        config["nanoplot"]["threads"]
    conda:
        "../envs/nanoplot.yml"
    shell:
        "NanoPlot "
            "--threads {threads} "
            "--alength " # aligned read length
            "--raw "
            "--bam {input} "
            "--outdir {output.outdir} "
            "> {log}"

# --------------------------------------------------------------------------- #
# Merge QC reports using MultiQC                                              #
# --------------------------------------------------------------------------- #
rule multiqc:
    """Merge QC reports using MultiQC"""
    message:
        "--- MultiQC"
    input:
        expand(rules.nanoplot_rawfastq.output.stats, sample = SAMPLES),
        expand(rules.nanoplot_filteredfastq.output.stats, sample = SAMPLES),
        expand(rules.nanoplot_aligned.output.stats, sample = SAMPLES)
    output:
        outdir = directory(os.path.join(outdir, "qc/multiqc")),
        report = os.path.join(outdir, "qc/multiqc/multiqc_report.html"),
        final = os.path.join(outdir, "variant_reports/MultiQC.html")
    log:
        os.path.join(outdir, "qc/multiqc/multiqc.log")
    threads:
        config["multiqc"]["threads"]
    params:
        qcdir = lambda x, output: os.path.dirname(str(output.outdir))
    conda:
        "../envs/multiqc.yml"
    shell:
        "multiqc "
            "--force "
            "--dirs "
            "--dirs-depth 2 "
            "--outdir {output.outdir} "
            "{params.qcdir} && "
        "cp {output.report} {output.final}"