# --------------------------------------------------------------------------- #
# Identification of structural variant using Sniffles                         #
# --------------------------------------------------------------------------- #
rule sniffles2:
    """Identification of structural variants using Sniffles2"""
    message:
        "--- SV calling using Sniffles2"
    input:
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
        bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
        reference=get_reference,
    output:
        os.path.join(outdir, "SV/sniffles2/{sample}.vcf"),
    log:
        os.path.join(outdir, "SV/sniffles2/logs/{sample}.log"),
    threads: config["sniffles2"]["threads"]
    conda:
        "../envs/sniffles2.yml"
    shell:
        "sniffles "
        "--threads {threads} "
        "--input {input.bam} "
        "--reference {input.reference} "
        "--vcf {output} "
        "{params} "
        "> {log}"


# --------------------------------------------------------------------------- #
# Identification of structural variant using cuteSV                           #
# --------------------------------------------------------------------------- #
rule cutesv:
    """Identification of structural variants using cuteSV"""
    message:
        "--- SV calling using cuteSV"
    input:
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
        bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
        reference=get_reference,
    output:
        vcfdir=directory(os.path.join(outdir, "SV/cutesv/{sample}")),
        vcf=os.path.join(outdir, "SV/cutesv/{sample}/{sample}.vcf"),
        finalvcf=os.path.join(outdir, "SV/cutesv/{sample}.vcf"),
    log:
        os.path.join(outdir, "SV/cutesv/logs/{sample}.log"),
    threads: config["cutesv"]["threads"]
    params:
        " ".join(config["cutesv"]["params"]),
    conda:
        "../envs/cutesv.yml"
    shell:
        "mkdir -p {output.vcfdir} && "
        "cuteSV "
        "{input.bam} "
        "{input.reference} "
        "{output.vcf} "
        "{output.vcfdir} "
        "--threads {threads} "
        "{params} "
        "2> {log} && "
        "cp {output.vcf} {output.finalvcf}"
