# --------------------------------------------------------------------------- #
# Identification of small nucleotide variants using Medaka                    #
# --------------------------------------------------------------------------- #
rule snv_medaka:
    """Identification of small nucleotide variants using Medaka"""
    message:
        "--- SNV calling using Medaka"
    input:
        reads = os.path.join(outdir, "filtered_reads/{sample}.fastq.gz"),
        reference = get_reference
    output:
        medakadir = directory(os.path.join(outdir, "SNV/medaka/{sample}")),
        vcf = os.path.join(outdir, "SNV/medaka/{sample}.vcf")
    log:
        stdout = os.path.join(outdir, "SNV/medaka/logs/{sample}.stdout"),
        stderr = os.path.join(outdir, "SNV/medaka/logs/{sample}.stderr")
    params:
        model = get_medaka_model
    threads:
        config["medaka"]["threads"]
    conda:
        "../envs/medaka.yml"
    shell:
        "medaka_haploid_variant "
            "-m {params.model} "
            "-t {threads} "
            "-o {output.medakadir} "
            "-i {input.reads} "
            "-r {input.reference} "
            "> {log.stdout} "
            "2> {log.stderr} && "
        "cp {output.medakadir}/medaka.annotated.vcf {output.vcf}"

# --------------------------------------------------------------------------- #
# Download models for Clair3                                                  #
# --------------------------------------------------------------------------- #
rule models_clair3:
    """Downloading basecalling models for Clair3"""
    message:
        "--- Downloading models for Clair3"
    output:
        directory(os.path.join(outdir, "SNV/clair3/model"))
    log:
        stdout = os.path.join(outdir, "SNV/clair3/logs/model.stdout"),
        stderr = os.path.join(outdir, "SNV/clair3/logs/model.stderr")
    params:
        download_model_for_clair3
    conda:
        "../envs/basic.yml"
    shell:
        "wget "
            "--directory-prefix {output} "
            "{params[0][0]} "
            "> {log.stdout} "
            "2> {log.stderr} && "
        "tar -xf {output}/{params[0][1]}.tar.gz -C {output} -v >> {log.stdout} && "
        "rm -rf {output}/{params[0][1]}.tar.gz && "
        "mv {output}/{params[0][1]}/* {output} && "
        "rm -rf {output}/{params[0][1]}"

# --------------------------------------------------------------------------- #
# Identification of small nucleotide variants using Clair3                    #
# --------------------------------------------------------------------------- #
rule snv_clair3:
    """Identification of small nucleotide variants using Clair3"""
    message:
        "--- SNV calling using Clair3"
    input:
        bam = os.path.join(outdir, "mapping/{sample}.bam"),
        bai = os.path.join(outdir, "mapping/{sample}.bam.bai"),
        reference = get_reference,
        model = rules.models_clair3.output
    output:
        clair3dir = directory(os.path.join(outdir, "SNV/clair3/{sample}")),
        vcf = os.path.join(outdir, "SNV/clair3/{sample}.vcf")
    log:
        stdout = os.path.join(outdir, "SNV/clair3/logs/{sample}.stdout"),
        stderr = os.path.join(outdir, "SNV/clair3/logs/{sample}.stderr")
    params:
        " ".join(config["clair3"]["params"])
    threads:
        config["clair3"]["threads"]
    conda:
        "../envs/clair3.yml"
    shell:
        "run_clair3.sh "
            "--threads {threads} "
            "{params} "
            "--bam_fn {input.bam} "
            "--ref_fn {input.reference} "
            "--model_path {input.model} "
            "--output {output.clair3dir} "
            "> {log.stdout} "
            "2> {log.stderr} && "
        "cp {output.clair3dir}/merge_output.vcf.gz {output.vcf}.gz && "
        "gunzip {output.vcf}.gz"