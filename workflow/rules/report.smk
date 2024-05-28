# --------------------------------------------------------------------------- #
# Collect VCF files for reporting                                             #
# --------------------------------------------------------------------------- #
rule collect_vcfs:
    """Collect region-filtered VCF files for reporting"""
    message:
        "--- Collect region-filtered VCF files for reporting"
    input:
        medaka=os.path.join(outdir, "SNV/medaka/{sample}.filtered.vcf"),
        clair3=os.path.join(outdir, "SNV/clair3/{sample}.filtered.vcf"),
        cutesv=os.path.join(outdir, "SV/cutesv/{sample}.filtered.vcf"),
        sniffles2=os.path.join(outdir, "SV/sniffles2/{sample}.filtered.vcf"),
    output:
        medaka=os.path.join(outdir, "variant_reports/{sample}/{sample}.medaka.vcf"),
        clair3=os.path.join(outdir, "variant_reports/{sample}/{sample}.clair3.vcf"),
        cutesv=os.path.join(outdir, "variant_reports/{sample}/{sample}.cutesv.vcf"),
        sniffles2=os.path.join(
            outdir, "variant_reports/{sample}/{sample}.sniffles2.vcf"
        ),
    log:
        os.path.join(outdir, "variant_reports/logs/{sample}.collect_vcfs.log"),
    conda:
        "../envs/basic.yml"
    shell:
        "cp {input.medaka} {output.medaka} > {log} && "
        "cp {input.clair3} {output.clair3} >> {log} && "
        "cp {input.cutesv} {output.cutesv} >> {log} && "
        "cp {input.sniffles2} {output.sniffles2} >> {log}"


# --------------------------------------------------------------------------- #
# Prepare VCF files for use with IGV                                          #
# --------------------------------------------------------------------------- #
rule prepare_vcfs:
    """Prepare VCF files for use with igv-reports"""
    message:
        "--- Preparing VCF files for IGV report"
    input:
        medaka=os.path.join(outdir, "variant_reports/{sample}/{sample}.medaka.vcf"),
        clair=os.path.join(outdir, "variant_reports/{sample}/{sample}.clair3.vcf"),
        cutesv=os.path.join(outdir, "variant_reports/{sample}/{sample}.cutesv.vcf"),
        sniffles=os.path.join(outdir, "variant_reports/{sample}/{sample}.sniffles2.vcf"),
    output:
        temp(os.path.join(outdir, "variant_reports/{sample}/all_variants.txt")),
    params:
        config["igv_region_length"],
    log:
        os.path.join(outdir, "variant_reports/logs/{sample}.prepare_vcfs.log"),
    conda:
        "../envs/basic.yml"
    script:
        "../scripts/merge_vcfs.py"


# --------------------------------------------------------------------------- #
# Generate IGV report                                                         #
# --------------------------------------------------------------------------- #
rule igv_reports:
    """Generate interactive IGV report"""
    message:
        "--- Generating IGV report"
    input:
        regions=os.path.join(outdir, "variant_reports/{sample}/all_variants.txt"),
        medaka=os.path.join(outdir, "variant_reports/{sample}/{sample}.medaka.vcf"),
        clair=os.path.join(outdir, "variant_reports/{sample}/{sample}.clair3.vcf"),
        cutesv=os.path.join(outdir, "variant_reports/{sample}/{sample}.cutesv.vcf"),
        sniffles=os.path.join(outdir, "variant_reports/{sample}/{sample}.sniffles2.vcf"),
        reference=get_reference,
        bam=os.path.join(outdir, "mapping/{sample}.bam"),
        bai=os.path.join(outdir, "mapping/{sample}.bam.bai"),
    output:
        os.path.join(outdir, "variant_reports/{sample}/{sample}_IGV.html"),
    log:
        os.path.join(outdir, "variant_reports/logs/{sample}.igv_reports.log"),
    params:
        gff=get_annotation,
    conda:
        "../envs/igv-reports.yml"
    shell:
        "create_report "
        "{input.regions} "
        "--sequence 1 "
        "--begin 2 "
        "--end 3 "
        "--fasta {input.reference} "
        "--tracks {params.gff} {input.bam} {input.medaka} {input.clair} {input.cutesv} {input.sniffles} "
        "--flanking 1000 "
        "--info-columns 'contig' 'variant region start' 'variant region end' 'variant details' "
        "--output {output} "
        "> {log}"


# --------------------------------------------------------------------------- #
# Generate report                                                             #
# --------------------------------------------------------------------------- #
rule report:
    """Generate report"""
    message:
        "--- Generating report"
    input:
        rawreadqc=os.path.join(outdir, "qc/raw_reads/{sample}/NanoPlot-data.tsv.gz"),
        filtreadqc=os.path.join(
            outdir, "qc/filtered_reads/{sample}/NanoPlot-data.tsv.gz"
        ),
        mappingstats=os.path.join(outdir, "mapping/{sample}.stats"),
        coverage=os.path.join(outdir, "mapping/{sample}.coverage"),
        five_plus=os.path.join(outdir, "mapping/{sample}.5ends.plus.counts"),
        five_minus=os.path.join(outdir, "mapping/{sample}.5ends.minus.counts"),
        three_plus=os.path.join(outdir, "mapping/{sample}.3ends.plus.counts"),
        three_minus=os.path.join(outdir, "mapping/{sample}.3ends.minus.counts"),
        medaka=os.path.join(outdir, "variant_reports/{sample}/{sample}.medaka.vcf"),
        clair3=os.path.join(outdir, "variant_reports/{sample}/{sample}.clair3.vcf"),
        cutesv=os.path.join(outdir, "variant_reports/{sample}/{sample}.cutesv.vcf"),
        sniffles=os.path.join(outdir, "variant_reports/{sample}/{sample}.sniffles2.vcf"),
        igv=os.path.join(outdir, "variant_reports/{sample}/{sample}_IGV.html"),
        multiqc=os.path.join(outdir, "variant_reports/MultiQC.html"),
    output:
        os.path.join(outdir, "variant_reports/{sample}/{sample}_overview.html"),
    log:
        os.path.join(outdir, "variant_reports/logs/{sample}.report.log"),
    params:
        maskedregions=get_bed_file_for_filtering,
        sharedvariants=get_shared_variants,
    conda:
        "../envs/report.yml"
    script:
        "../scripts/report.Rmd"
