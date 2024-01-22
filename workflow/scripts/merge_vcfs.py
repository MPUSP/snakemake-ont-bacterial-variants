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

def determine_svtype(short):
    types = {
        "DEL": "deletion",
        "INS": "insertion",
        "BND": "translocation breakend",
        "DUP": "duplication",
        "INV": "inversion"}
    return(types[short])

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

# Select regions for resulting IGV report
initiate = True
regions = []
for index, row in df.iterrows():
    currentcontig = row["CHROM"]
    currentposition = row["POS"]
    if row["ID"].startswith("Sniffles"):
        vartype = determine_svtype(row["ID"].split(".")[1])
    elif row["ID"].startswith("cuteSV"):
        vartype = determine_svtype(row["ID"].split(".")[1])
    else:
        if len(row["ALT"]) > len(row["REF"]): vartype = "insertion"
        elif len(row["ALT"]) < len(row["REF"]): vartype = "deletion"
        else: vartype = "SNV"
    currentinfo = F"{vartype} at pos. {row['POS']} in {row['CHROM']} called by {row['TOOL']}"
    if initiate:
        entry = [currentcontig, currentposition, currentposition, currentinfo]
        previouscontig = currentcontig
        previousposition = currentposition
        initiate = False
    else:
        currentcontig = row["CHROM"]
        currentposition = row["POS"]
        if (currentcontig == previouscontig) and (int(currentposition) - int(previousposition) <= 1000):
            entry[2] = currentposition
            entry[3] = entry[3] + "; " + currentinfo
        else:
            regions.append(entry)
            entry = [currentcontig, currentposition, currentposition, currentinfo]
            previouscontig = currentcontig
            previousposition = currentposition
regions.append(entry)

# Write to outfile
regions = pd.DataFrame(regions, columns = ["contig", "variant region start", "variant region end", "variant details"])
regions.to_csv(str(snakemake.output), sep = "\t", index = False)