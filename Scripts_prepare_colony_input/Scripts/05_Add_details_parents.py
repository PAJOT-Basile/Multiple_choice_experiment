#!/bin/python3

###############################
########## Libraries ##########
###############################
import pandas as pd
import yaml
import pysam as vcf
print("\t\tPreparing Parents' stats")

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


def potential_parents(samples, metadata, sex):
    return ([ID for ID in metadata["ID_DNA_RAD"][(
        metadata["Family_level"] == "parents") & (metadata["Sex"] == sex)] if ID in samples])


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

# Import metadata
metadata = pd.read_csv(config_file["Metadata"], sep="\t", header=0)

###############################
####### Compute stats #########
###############################
# Get the samples from the vcf
samples_from_vcf = get_samples_from_vcf(vcf_file)

# Get the names of the potential parents
potential_fathers = potential_parents(samples_from_vcf, metadata, "M")
potential_mothers = potential_parents(samples_from_vcf, metadata, "F")

# Get the numbers of potential fathers and mothers
nb_fathers, nb_mothers = len(potential_fathers), len(potential_mothers)

# Get the number of broods
nb_broods = len(
    metadata["Mother_ID"][metadata["Family_level"] == "offspring"].value_counts())

# Compute the probability of the dad and mom be included in the candidates. This is calculated as follows:
# number of broods divided by the number of potential dads if we consider that each male contributed to at least on brood, wich is not true.
# number of potential females divided by the number of broods if we consider that all females produced some offspring
prop_fathers, prop_mothers = nb_broods/nb_fathers, nb_mothers/nb_broods

with open(outfile, "a") as f:
    f.write(f"""{round(prop_fathers, 1)}  {round(prop_mothers, 1)} \t !prob. of dad/mum included in the candidates
{nb_fathers}  {nb_mothers} \t !numbers of candidate males and females
""")
