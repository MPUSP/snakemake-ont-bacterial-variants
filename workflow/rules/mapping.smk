# --------------------------------------------------------------------------- #
# Read filtering                                                              #
# --------------------------------------------------------------------------- #
rule filtlong:
    """Filtering of nanopore reads"""
    message:
        "--- Read filtering"
    input:
        get_fastq
    output:
        os.path.join(outdir, "filtered_reads/{sample}.fastq.gz")
    log:
        os.path.join(outdir, "filtered_reads/logs/{sample}.log")
    params:
        default = " ".join(config["filtlong"]["params"])
    conda:
        "../envs/filtlong.yml"
    shell:
        "filtlong "
            "{params.default} "
            "{input} "
            "2> {log} | "
        "gzip > {output}"

# --------------------------------------------------------------------------- #
# Mapping against reference genome                                            #
# --------------------------------------------------------------------------- #
rule mapping:
    """Mapping of ONT data against reference genome"""
    message:
        "--- Mapping against reference genome"
    input:
        reads = os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
        reference = get_reference
    output:
        sam = temp(os.path.join(outdir, "mapping/{sample}.sam")),
        bam = os.path.join(outdir, "mapping/{sample}.bam"),
        bai = os.path.join(outdir, "mapping/{sample}.bam.bai"),
        stats = temp(os.path.join(outdir, "mapping/{sample}.stats"))
    log:
        ngmlr = os.path.join(outdir, "mapping/logs/{sample}.ngmlr.log"),
        samtools = os.path.join(outdir, "mapping/logs/{sample}.samtools.log")
    params:
        tmpdir = os.path.join(outdir, "mapping/{sample}")
    threads:
        config["ngmlr"]["threads"]
    conda:
        "../envs/ngmlr.yml"
    shell:
        "ngmlr "
            "-t {threads} "
            "-x ont "
            "-r {input.reference} "
            "-q {input.reads} "
            "-o {output.sam} "
            "2> {log.ngmlr} && "
        "samtools sort "
            "{output.sam} "
            "-@ {threads} "
            "-T {params.tmpdir} "
            "-O bam "
            "> {output.bam} "
            "2> {log.samtools} && "
        "samtools index {output.bam} 2>> {log.samtools} && "
        "samtools stats {output.bam} > {output.stats} 2>> {log.samtools}"

# --------------------------------------------------------------------------- #
# Genome coverage                                                             #
# --------------------------------------------------------------------------- #
rule genomecoverage:
    """Generation of nucleotide-precise genome coverage depths"""
    message:
        "--- Genome coverage"
    input:
        os.path.join(outdir, "mapping/{sample}.bam")
    output:
        temp(os.path.join(outdir, "mapping/{sample}.coverage"))
    log:
        os.path.join(outdir, "mapping/logs/{sample}.genomecoverage.log")
    conda:
        "../envs/ngmlr.yml"
    shell:
        "bedtools genomecov "
            "-d "
            "-ibam {input} "
            "> {output} "
            "2> {log}"

# --------------------------------------------------------------------------- #
# Alignment ends                                                              #
# --------------------------------------------------------------------------- #
rule alignmentends:
    """Generation of 5' and 3' alignment end distributions by strand"""
    message:
        "--- Alignment end distribution"
    input:
        os.path.join(outdir, "mapping/{sample}.bam")
    output:
        five_plus = temp(os.path.join(outdir, "mapping/{sample}.5ends.plus.counts")),
        five_minus = temp(os.path.join(outdir, "mapping/{sample}.5ends.minus.counts")),
        three_plus = temp(os.path.join(outdir, "mapping/{sample}.3ends.plus.counts")),
        three_minus = temp(os.path.join(outdir, "mapping/{sample}.3ends.minus.counts"))
    log:
        os.path.join(outdir, "mapping/logs/{sample}.alignmentends.log")
    conda:
        "../envs/ngmlr.yml"
    shell:
        "bedtools genomecov -d -ibam {input} "
            "-5 -strand + > {output.five_plus} 2> {log} && "
        "bedtools genomecov -d -ibam {input} "
            "-5 -strand - > {output.five_minus} 2>> {log} && "
        "bedtools genomecov -d -ibam {input} "
            "-3 -strand + > {output.three_plus} 2>> {log} && "
        "bedtools genomecov -d -ibam {input} "
            "-3 -strand - > {output.three_minus} 2>> {log}"

