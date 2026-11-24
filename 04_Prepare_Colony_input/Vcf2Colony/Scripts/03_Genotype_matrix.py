#!/bin/python3


###############################
########## Libraries ##########
###############################
from warnings import simplefilter
import string
import random as rd
import numpy as np
import pandas as pd
import yaml
import pysam as vcf
from math import ceil
import argparse
simplefilter(action="ignore", category=pd.errors.PerformanceWarning)

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


def sample_list_with_list_indices(list_to_sample, list_of_indices):
    # ' Sample a list using a list with the indices to sample
    return ([list_to_sample[sample_index]
            for sample_index in list_of_indices])


def randomly_chose_genotypes():
    # ' Sample random genotypes in the vcf format (either 0 or 2)
    random_genotype = round(rd.uniform(0, 1)) * 2
    return ([str(random_genotype), str(2 - random_genotype)])


def standardise_for_colony(x):
    # ' Standardise a value from 0-1 to 1-2 and add 0 for missing data
    try:
        return int(x) + 1
    except ValueError:
        # If there is an error, we are missing a value here, so the default is 0
        return (0)


def prepare_genotype_for_colony(genotype):
    # Then transform the genotypes for colony use (standardise if between 1 and 2 compared to vcf which is between 0 and 2)
    diploid_genotypes = [standardise_for_colony(
        haploid_genotype) for haploid_genotype in genotype.split("/")]
    return diploid_genotypes


def transform_genotypes_for_colony(genotypes, sample_indices):
    # ' Transform the vcf genotypes into colony-friendly genotypes
    # First, sample the individuals that we want
    genotypes_of_samples = sample_list_with_list_indices(
        genotypes, sample_indices)
    # Then standardise the genotypes for colony
    genotypes_colony = [prepare_genotype_for_colony(
        genotype) for genotype in genotypes_of_samples]
    return (genotypes_colony)


def get_indices_samples_vcf(vcf_file, samples):
    i = 1
    # Read only the first line of the vcf file
    for rec in vcf_file.fetch():
        if i > 1:
            break
        i += 1
    # Return the indices of the samples
    return ([list(rec.samples).index(i)
            for i in rec.samples if i in samples])


def prepare_write_in_outfile(genotypes, id):
    genotypes_to_write = "\t".join([str(genotype[id])
                                   for genotype in genotypes])
    return (genotypes_to_write)


###############################
###### Import variables #######
###############################
# Import command line arguments
parser = argparse.ArgumentParser(
    prog="Get_genotype_matrix",
    description="Gets the genotype matrix for a sub-sample of individuals"
)
parser.add_argument("-s", "--sub_sample", help="Add a filtration condition on the "
                    "metadata. The default is no filtration and the full genotype matrix will be returned. "
                    "The possible filtrations are: 'Offspring' to select only offsprings, 'Parents_male' to "
                    "select only male parents and 'Parents_female' to select only female parents.")
args = parser.parse_args()
sub_sample = string.capwords(args.sub_sample)

print("".join(["\tMaking genotype matrix for ", sub_sample]))
# Import variables in the configuration file
with open("Data/Input.yaml", "r") as info:
    config_file = yaml.safe_load(info)

# Get the name of the output file
temp_file = "".join([config_file["Tmp_dir"], sub_sample, "_genotypes.txt"])
rd.seed(config_file["Random_seed"])
###############################
######## Import data ##########
###############################
# Import the vcf file
vcf_file = vcf.VariantFile(config_file["Input_file"])

# Import metadata
metadata = pd.read_csv(config_file["Metadata"], sep="\t", header=0)

# Check if there is some filtration to do on the metadata using the imported parsed arguments
match sub_sample:
    case "Offspring" | "Offsprings":
        metadata = metadata[(metadata["Family_level"] == "offspring")]
    case "Parents_male" | "Parents_males" | "Parents_Males" | "Parents_Male":
        metadata = metadata[(metadata["Family_level"] ==
                             "parents") & (metadata["Sex"] == "M")]
    case "Parents_female" | "Parents_females" | "Parents_Females" | "Parents_Female":
        metadata = metadata[(metadata["Family_level"] ==
                             "parents") & (metadata["Sex"] == "F")]
    case _:
        pass

###############################
###### Genotype matrix ########
###############################
# Make a list of all the samples to keep
samples = [sample for sample in get_samples_from_vcf(
    vcf_file) if sum(metadata["ID_DNA_RAD"].str.contains(sample)) > 0]
# Open the output file to write line by line
tmp_out_file = open(temp_file, "w+")
# Add the name of the samples
tmp_out_file.write("\t".join(samples) + "\n")

# Get the indices of the samples in the samples of the vcf
sample_indices = get_indices_samples_vcf(vcf_file, samples)

# Iterate over the lines of the vcf and print them to the file
for rec in vcf_file.fetch():

    # Get the genotypes for all individuals from the vcf line
    genotypes = [information.split(":")[0]
                 for information in str(rec).split("\t") if "/" in information]
    # Sample genotypes for the indivdiuals we want and prepare it for colony
    # (standardize values of genotypes between 1 and 2 against 0 and 2 for vcf and
    # add 0 for missing values)
    sample_genotypes_for_colony = np.array(transform_genotypes_for_colony(
        genotypes, sample_indices))
    # Prepare the array of genotypes for them to be written in the output file
    genotypes_to_write = "\n".join([
        prepare_write_in_outfile(sample_genotypes_for_colony, 0),
        prepare_write_in_outfile(sample_genotypes_for_colony, 1)
    ]) + "\n"
    # Write in the output file
    tmp_out_file.write(genotypes_to_write)

# Close the connection to the output file
tmp_out_file.close()
