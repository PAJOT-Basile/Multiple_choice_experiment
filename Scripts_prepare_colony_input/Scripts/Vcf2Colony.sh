#!/bin/bash

#SBATCH --partition=fast
#SBATCH --account=sexisol
#SBATCH --job-name=Vcf2Colony
#SBATCH --mem=50G
#SBATCH --cpus-per-task=20

# Create the conda environment
# source Scripts/00_Load_conda_env.sh

module load conda
# Activate the created conda environment
source activate ./Data/envs/Vcf2Colony/

# Start running the scripts

# Import some variables to be used here
output_file=$(grep "Outfolder" Data/Input.yaml | cut -d" " -f2 | perl -pe 's/"//g')colony.dat
temp_folder=$(grep "Outfolder" Data/Input.yaml | cut -d" " -f2 | perl -pe 's/"//g')
# Create temporary directory if it does not exist
mkdir -p "${temp_folder}"


python3 Scripts/01_Prepare_header.py
Rscript Scripts/02_Table_marker_type.r
echo -e "\n" >> $output_file

for samples in $(echo "Offspring Parents_male Parents_female"); do
   python3 Scripts/03_Genotype_matrix.py -s "${samples}"
   source Scripts/04_Transpose_file.sh -s "${samples}"
   if [[ "${samples}" == "Offspring" ]]; then
      python3 Scripts/05_Add_details_parents.py
      echo -e "\n" >> "${output_file}"
   fi
done


echo "Done!"
conda deactivate
module unload conda