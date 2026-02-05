#!/bin/python3

###############################
########## Libraries ##########
###############################
import pandas as pd
import yaml
import pysam as vcf
print("Preparing parent-offspring sibship tables")

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


def get_number_offspring_mother_dyads(metadata):
    list_mother_ids = metadata["Label"][(
        metadata["Family_level"] == "parents") & (metadata["Sex"] == "F")].values
    metadata_offspring = metadata[(metadata["Family_level"] == "offspring")]
    offspring_to_count_in_dyad = [row["ID_DNA_RAD"] for _, row in metadata_offspring.iterrows(
    ) if row["Mother_ID"] in list_mother_ids]
    return (len(offspring_to_count_in_dyad))


def get_mother_offspring_dyads(mother_ID, metadata):
    mothers_offspring = metadata[(metadata["Mother_ID"] == mother_ID)]
    table_mother_rad = metadata["ID_DNA_RAD"][(
        metadata["Label"] == mother_ID)].values
    if table_mother_rad.shape[0] == 0:
        return None
    mother_rad_id = str(table_mother_rad[0])
    mother_offspring_dyads = pd.DataFrame({
        "Offspring": mothers_offspring["ID_DNA_RAD"],
        "Mother_ID": mother_rad_id
    })
    return (mother_offspring_dyads)


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

# TODO: test the values of the paternity/maternity exclusion threshold (start at 0):
# def = provide the number of offspring in the offspring samples that have known father sincluded in the
# candidate male samples and the maximum #mismatches allowed for a known father-offspring dyad.
Exclusion_threshold = config_file["Paternity_and_maternity_exclusion_threshold"]

################################################
######## Get the mother-offspring dyads ########
################################################
# Get the number of mother-offsrping dyads
Nb_mother_offspring_dyads = get_number_offspring_mother_dyads(
    metadata=metadata)

# Get the list of all the mothers
all_mothers = metadata.Mother_ID.dropna().unique()


###############################
######## Write output #########
###############################
with open(outfile, "a") as out:
    out.write(
        f"""0  {Exclusion_threshold} \t !#known father-offspring dyads, paternity exclusion threshold
\n
{Nb_mother_offspring_dyads}  {Exclusion_threshold} \t !#known mother-offspring dyads, maternity exclusion threshold
"""
    )

# Open the output file to write the dyads inside
out_file = open(outfile, "a")

# Iterate over the mothers to get the mother-offspring dyads
for mother in all_mothers:
    mother_offspring_dyads = get_mother_offspring_dyads(
        mother_ID=mother, metadata=metadata)
    if mother_offspring_dyads is None:
        continue
    for _, row in mother_offspring_dyads.iterrows():
        out_file.write("\t".join(row.values) + "\n")

# Close the output file
out_file.close()

# Add the rest of the known sibships in the table
with open(outfile, "a") as out:
    out.write(
        f"""\n
0 \t !#known paternal sibship with unknown fathers
\n
0 \t !#known maternal sibship with unknown mothers
\n
0 \t !#known paternity exclusions
\n
0 \t !#known maternity exclusions
\n
0 \t !#known paternal sibship exclusions
\n
0 \t !#known maternal sibship exclusions
"""
    )
