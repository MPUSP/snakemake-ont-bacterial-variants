import os
import pandas as pd

# --------------------------------------------------------------------------- #
# General helper functions                                                    #
# --------------------------------------------------------------------------- #
def read_sample_sheet():
    samples = pd.read_csv(
        config["samples"],
        sep="\t"
    ).set_index("sample", drop = False)
    return(samples)

def list_reference_genomes():
    paths = SAMPLEINFO["reference"].unique()
    genomes = {}
    for path in paths:
        ident = ".".join(path.split('/')[-1].split(".")[:-1])
        if not ident in genomes.keys():
            genomes[ident] = path
    return(genomes)

def generate_results(wildcards):
    return([
        os.path.join(outdir, F"variant_reports/{sample}/{sample}_overview.html")
        for sample in SAMPLES
    ])

# --------------------------------------------------------------------------- #
# Helper functions using wildcards                                            #
# --------------------------------------------------------------------------- #
def get_fastq(wildcards):
    return(SAMPLEINFO.loc[wildcards.sample, "fastq"])

def get_reference(wildcards):
    return(SAMPLEINFO.loc[wildcards.sample, "reference"])

def get_annotation(wildcards):
    if not pd.isna(SAMPLEINFO.loc[wildcards.sample, "annotation"]):
        if os.path.exists(SAMPLEINFO.loc[wildcards.sample, "annotation"]):
            return(SAMPLEINFO.loc[wildcards.sample, "annotation"])
    else:
        return("")

def download_model_for_clair3(wildcards):
    path2model = {
        "r1041_e82_400bps_sup_v4.2.0": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v420.tar.gz",
        "r1041_e82_400bps_sup_v4.3.0": "https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/r1041_e82_400bps_sup_v430.tar.gz"
    }
    return(
        path2model[config["basecalling_model"]],
        path2model[config["basecalling_model"]].split('/')[-1].split(".tar.gz")[0])

def get_medaka_model(wildcards):
    modeldetails = config["basecalling_model"].split("_")
    modeldetails.insert(-1, "variant")
    return("_".join(modeldetails))

def get_vcf(wildcards):
    assigntype = {"medaka": "SNV", "clair3": "SNV", "cutesv": "SV", "sniffles2": "SV"}
    return(os.path.join(outdir, assigntype[wildcards.tool], "{tool}/{sample}.vcf"))

def get_bed_file_for_filtering(wildcards):
    return(str(SAMPLEINFO.loc[wildcards.sample, "masked_regions"]))

def get_shared_variants(wildcards):
    if config["remove_common_variants"] and len(SAMPLES) > 1:
        return([
            str(os.path.join(outdir, "SNV/clair3/common_variants.vcf")),
            str(os.path.join(outdir, "SNV/medaka/common_variants.vcf")),
            str(os.path.join(outdir, "SV/cutesv/common_variants.vcf")),
            str(os.path.join(outdir, "SV/sniffles2/common_variants.vcf"))
        ])
    else:
        return()