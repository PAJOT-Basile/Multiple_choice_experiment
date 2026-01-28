#!/bin/python3
# TODO: the number of offspring in the sample and the mean number of offspring per mother and father do not use the vcf input file, to tune
###############################
########## Libraries ##########
###############################
import pandas as pd
import yaml
import pysam as vcf
print("Preparing header")

###############################
###### Import variables #######
###############################
with open("Data/Input.yaml", "r") as info:
    config_file = yaml.safe_load(info)


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
# Get the number of offsrping
nb_offspring = len(metadata.query("Family_level == 'offspring'").index)

# Get the mean number of offspring per female
nb_offspring_per_female = round(metadata.groupby("Mother_ID").size().mean(), 2)

# Get the number of loci
nb_loci = 0
nb_alleles_per_loc = []
for rec in vcf_file.fetch():
    nb_loci += 1
    nb_alleles_per_loc.append(len(rec.alleles))


###############################
#### Prepare output header ####
###############################
with open("".join([config_file["Outfolder"], "colony.dat"]), "w") as f:
    f.write(f"""'{config_file["Output_name"]}' \t !Dataset name
'{config_file["Output_name"]}' \t !Output file name
{nb_offspring} \t ! Number of offspring in the sample
{nb_loci} \t ! Number of loci
{config_file["Random_seed"]}                   ! Seed for random number generator
{config_file["Update_allele_freq"]} \t ! 0/1=Not updating/updating allele frequency
{config_file["Dioecious_or_monoecious"]} \t ! 2/1=Dioecious/Monoecious species
{config_file["Inbreeding"]} \t ! 0/1=No inbreeding/inbreeding
{config_file["Ploidy"]} \t ! 0/1=Diploid species/HaploDiploid species
{config_file["Gamy"][0]}  {config_file["Gamy"][1]} \t ! 0/1=Polygamy/Monogamy for males & females
{config_file["Clone"]} \t ! 0/1=Clone inference =No/Yes
{config_file["Sibship_size_scaling"]} \t ! 0/1=Full sibship size scaling =No/Yes
{config_file["Sibship_size_prior"]} {nb_offspring_per_female} {nb_offspring_per_female} \t ! 0,1,2,3=No,weak,medium,strong sibship size prior; mean paternal & meteral sibship size
{config_file["Known_allele_freq"]} \t ! 0/1=Unknown/Known population allele frequency
{" ".join([str(x) for x in nb_alleles_per_loc])} \t !Number of alleles per locus
{config_file["Nb_runs"]} \t ! Number of runs
{config_file["Length_run"]} \t ! Length of run
{config_file["Monitor_method_by_iterate_number_in_second"]} \t ! 0/1=Monitor method by Iterate#/Time in second
{config_file["Monitor_interval_in_iterate_number_in_second"]} \t ! Monitor interval in Iterate# / in seconds
{config_file["Windows_verion"]} \t ! non-Windows version
{config_file["Fulllikelihood"]} \t ! Fulllikelihood
{config_file["Precision_for_likelihood"]} \t ! 1/2/3=low/medium/high Precision for Fulllikelihood

""")
