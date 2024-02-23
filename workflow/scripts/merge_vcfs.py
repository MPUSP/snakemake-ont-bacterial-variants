#!/usr/bin/env python
import pandas as pd

def read_vcf(vcffile):
    out = []
    columnnames = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                   "INFO", "FORMAT", "SAMPLE"]
    with open(vcffile, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line.startswith("#"):
                out.append(line.split("\t"))
    return(pd.DataFrame(out, columns = columnnames))

def determine_type(series):
    types = {
        "DEL": "deletion",
        "INS": "insertion",
        "BND": "translocation breakend",
        "DUP": "duplication",
        "INV": "inversion"}
    if series["ID"].startswith("Sniffles"):
        vartype = types[series["ID"].split(".")[1]]
    elif series["ID"].startswith("cuteSV"):
        vartype = types[series["ID"].split(".")[1]]
    else:
        if len(series["ALT"]) > len(series["REF"]): vartype = "insertion"
        elif len(series["ALT"]) < len(series["REF"]): vartype = "deletion"
        else: vartype = "SNV"
    return(vartype)

# Import VCF files
medaka = read_vcf(snakemake.input["medaka"])
medaka["TOOL"] = ["Medaka" for i in range(medaka.shape[0])]
clair = read_vcf(snakemake.input["clair"])
clair["TOOL"] = ["Clair3" for i in range(clair.shape[0])]
sniffles = read_vcf(snakemake.input["sniffles"])
sniffles["TOOL"] = ["Sniffles" for i in range(sniffles.shape[0])]
cutesv = read_vcf(snakemake.input["cutesv"])
cutesv["TOOL"] = ["cuteSV" for i in range(cutesv.shape[0])]

# Merge and sort variant information
df = pd.concat([medaka, clair, sniffles, cutesv], ignore_index = True)
df["POS"] = pd.to_numeric(df["POS"])
df = df.sort_values(by = ["CHROM", "POS"])
df["TYPE"] = df.apply(determine_type, axis = 1)

# Select regions for resulting IGV report
distance_threshold = int(snakemake.params)
initiate = True
regions = []
for contig in df["CHROM"].unique():
    currentcontig = contig
    for position in df[df["CHROM"] == contig]["POS"].unique():
        for vartype in df[(df["CHROM"] == contig) & (df["POS"] == position)]["TYPE"].unique():
            subset = df[(df["CHROM"] == contig) & (df["POS"] == position) & (df["TYPE"] == vartype)]
            temp = []
            for row in subset.index:
                temp.append(F"{subset.loc[row, 'TOOL']}: quality of {subset.loc[row, 'QUAL']}")
            vardetails = F"{vartype} at pos. {position} ({'; '.join(temp)})"
        currentposition = position
        if initiate:
            entry = [currentcontig, currentposition, currentposition, vardetails]
            previouscontig = currentcontig
            previousposition = currentposition
            initiate = False
        if (currentcontig == previouscontig) and (currentposition - previousposition) <= distance_threshold:
                entry[2] = currentposition
                entry[3] = entry[3] + "; " + vardetails
        else:
            regions.append(entry)
            entry = [currentcontig, currentposition, currentposition, vardetails]
            previouscontig = currentcontig
            previousposition = currentposition
regions.append(entry)

# Write to outfile
regions = pd.DataFrame(regions, columns = ["contig", "variant region start", "variant region end", "variant details"])
regions.to_csv(str(snakemake.output), sep = "\t", index = False)