#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=sexisol
#SBATCH --job-name=Vcf2Colony
#SBATCH --mem=10G
#SBATCH --cpus-per-task=1

# module load conda
# First, check if the conda environment is created. If not, create it
if [[ ! -d "Data/envs/Vcf2Colony" ]]; then
   conda env create -y --file Data/conda_env.yaml -p Data/envs/Vcf2Colony
   conda config -p ./Data/envs/Vcf2Colony
   conda init --all --user
else
   echo "Vcf2Colony conda environment already exists"
fi


# module unload conda