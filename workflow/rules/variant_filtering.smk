if config["remove_common_variants"] and len(SAMPLES) > 1:

    # --------------------------------------------------------------------------- #
    # Compression and indexing of VCF files for use with bcftools                 #
    # --------------------------------------------------------------------------- #
    rule prepare4bcftools:
        """Compress and index all VCF files"""
        message:
            "--- Identification of common variants: VCF file indexing"
        input:
            os.path.join(outdir, "{vartype}/{tool}/{sample}.vcf"),
        output:
            compressed=temp(os.path.join(outdir, "{vartype}/{tool}/{sample}.vcf.gz")),
            indexed=temp(os.path.join(outdir, "{vartype}/{tool}/{sample}.vcf.gz.tbi")),
        log:
            stdout=os.path.join(
                outdir, "{vartype}/{tool}/logs/{sample}_indexing.stdout"
            ),
            stderr=os.path.join(
                outdir, "{vartype}/{tool}/logs/{sample}_indexing.stderr"
            ),
        conda:
            "../envs/bcftools.yml"
        shell:
            "bcftools sort {input} -Oz -o {output.compressed} > {log.stdout} 2> {log.stderr} && "
            "bcftools index -f --tbi {output.compressed} -o {output.indexed} >> {log.stdout} 2>> {log.stderr}"

    # --------------------------------------------------------------------------- #
    # Identification of common shared variants across samples                     #
    # --------------------------------------------------------------------------- #
    rule bcftools:
        """Identification of common variants shared by all samples"""
        message:
            "--- Identification of common variants: variant detection"
        input:
            compressed=expand(
                os.path.join(outdir, "{{vartype}}/{{tool}}/{sample}.vcf.gz"),
                sample=SAMPLES,
            ),
            indexed=expand(
                os.path.join(outdir, "{{vartype}}/{{tool}}/{sample}.vcf.gz.tbi"),
                sample=SAMPLES,
            ),
        output:
            os.path.join(outdir, "{vartype}/{tool}/common_variants.vcf"),
        log:
            stdout=os.path.join(outdir, "{vartype}/{tool}/logs/bcftools.stdout"),
            stderr=os.path.join(outdir, "{vartype}/{tool}/logs/bcftools.stderr"),
        params:
            vcf_count=lambda x, input: len([str(elem) for elem in input.compressed]),
        conda:
            "../envs/bcftools.yml"
        shell:
            "bcftools isec "
            "-n={params.vcf_count} "
            "-w1 "
            "{input.compressed} "
            "-o {output}"

    # --------------------------------------------------------------------------- #
    # Convert VCF file with common variants to TXT                                #
    # --------------------------------------------------------------------------- #
    rule preparevcf4vcftools:
        """Conversion of VCF file with common variants to tab-delimited list"""
        message:
            "--- Identification of common variants: VCF-to-TXT conversion"
        input:
            os.path.join(outdir, "{vartype}/{tool}/common_variants.vcf"),
        output:
            os.path.join(outdir, "{vartype}/{tool}/common_variants.txt"),
        shell:
            """awk '!/^#/ {{print $1 "\t" $2}}' {input} > {output}"""

    # --------------------------------------------------------------------------- #
    # Filter variants                                                             #
    # --------------------------------------------------------------------------- #
    rule filter_variants_with_common:
        """Filtering of variants"""
        message:
            "--- Filtering variants"
        input:
            vcf=os.path.join(outdir, "{vartype}/{tool}/{sample}.vcf"),
            commonvars=os.path.join(outdir, "{vartype}/{tool}/common_variants.txt"),
        output:
            os.path.join(outdir, "{vartype}/{tool}/{sample}.filtered.vcf"),
        log:
            os.path.join(outdir, "{vartype}/{tool}/logs/{sample}_filtering.log"),
        params:
            bed=get_bed_file_for_filtering,
            prefix=lambda x, input: os.path.splitext(str(input.vcf))[0],
            minq=lambda wildcards: config["quality_threshold"][wildcards.tool],
        conda:
            "../envs/vcftools.yml"
        shell:
            "if [ -f '{params.bed}' ]; then "
            "vcftools "
            "--exclude-positions {input.commonvars} "
            "--minQ {params.minq} "
            "--vcf {input.vcf} "
            "--recode "
            "--keep-INFO-all "
            "--stdout "
            "2> {log} | "
            "vcftools "
            "--exclude-bed {params.bed} "
            "--vcf - "
            "--recode "
            "--keep-INFO-all "
            "--out {params.prefix} "
            "2>> {log} && "
            "mv {params.prefix}.recode.vcf {output}; "
            "else "
            "vcftools "
            "--exclude-positions {input.commonvars} "
            "--minQ {params.minq} "
            "--vcf {input.vcf} "
            "--recode "
            "--keep-INFO-all "
            "--out {params.prefix} "
            "2> {log} && "
            "mv {params.prefix}.recode.vcf {output} && "
            "echo 'No BED file for filtering provided' > {log}; "
            "fi "

else:

    # --------------------------------------------------------------------------- #
    # Filter variants                                                             #
    # --------------------------------------------------------------------------- #
    rule filter_variants_without_common:
        """Filtering of variants"""
        message:
            "--- Filtering variants"
        input:
            vcf=os.path.join(outdir, "{vartype}/{tool}/{sample}.vcf"),
        output:
            os.path.join(outdir, "{vartype}/{tool}/{sample}.filtered.vcf"),
        log:
            os.path.join(outdir, "{vartype}/{tool}/logs/{sample}_filtering.log"),
        params:
            bed=get_bed_file_for_filtering,
            prefix=lambda x, input: os.path.splitext(str(input.vcf))[0],
            minq=lambda wildcards: config["quality_threshold"][wildcards.tool],
        conda:
            "../envs/vcftools.yml"
        shell:
            "if [ -f '{params.bed}' ]; then "
            "vcftools "
            "--exclude-bed {params.bed} "
            "--minQ {params.minq} "
            "--vcf {input.vcf} "
            "--recode "
            "--keep-INFO-all "
            "--out {params.prefix} "
            "2> {log} && "
            "mv {params.prefix}.recode.vcf {output}; "
            "else "
            "vcftools "
            "--minQ {params.minq} "
            "--vcf {input.vcf} "
            "--recode "
            "--keep-INFO-all "
            "--out {params.prefix} "
            "2> {log} && "
            "mv {params.prefix}.recode.vcf {output} && "
            "echo 'No BED file for filtering provided' > {log}; "
            "fi "
