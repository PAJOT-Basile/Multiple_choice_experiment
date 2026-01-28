#!/bin/python3

###############################
########## Libraries ##########
###############################
import pandas as pd
import yaml
import pysam as vcf
print("Preparing Parents' stats")

###############################
###### Useful functions #######
###############################


def get_samples_from_vcf(vcf_file):
    i = 1
    for rec in vcf_file.fetch():
        if i > 1:
            break
        samples = rec.samples
        i += 1
    return (list(samples))


def filter_metadata_on_vcf_samples(samples, metadata):
    return ([ID for ID in metadata["ID_DNA_RAD"] if ID in samples])


###############################
###### Import variables #######
###############################
with open("Data/Input.yaml", "r") as info:
    config_file = yaml.safe_load(info)

outfile = "".join([config_file["Outfolder"], "colony.dat"])
###############################
######## Import data ##########
###############################
# Import the vcf file
vcf_file = vcf.VariantFile(config_file["Input_file"])

sequenced_samples = get_samples_from_vcf(vcf_file)

# Import metadata
metadata = pd.read_csv(config_file["Metadata"], sep="\t", header=0)
# Filter the metadata to keep only the individuals that were sequenced
metadata = metadata[metadata["ID_DNA_RAD"].isin(sequenced_samples)]


###############################
######## Write output #########
###############################
with open(outfile, "a") as out:
    out.write(
        f"""0  0 \t !#known fater-offspring dyads, paternity exclusion threshold
0  0 \t !#known moter-offspring dyads, maternity exclusion threshold
"""
    )
