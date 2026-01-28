#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=sexisol
#SBATCH --job-name=Vcf2Colony
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1

# Import the directories in which to get and save the files
output_file=$(grep "Outfolder" Data/Input.yaml | cut -d" " -f2 | perl -pe 's/"//g')colony.dat
temp_folder=$(grep "Tmp_dir" Data/Input.yaml | cut -d" " -f2 | perl -pe 's/"//g')

# Parse the arguments of the file
while getopts s: opt; do
    case "${opt}" in
    s)
        SAMPLES="${OPTARG}"
        ;;
    esac
done

# Get the name of the input file
infile="${temp_folder}${SAMPLES}_genotypes.txt"
# Count the number of columns
cols=`head -n1 "${infile}" | wc -w`
# Iterate over columns and transform them into lines
for (( i=1; i<="${cols}"; i++)); do
   cut -f $i "${infile}" | tr $'\n' $'\t' | sed -e "s/\t$/\n/g" | perl -pe "s/^/\t/g">> "${output_file}"
done

echo -e "\n" >> "${output_file}"