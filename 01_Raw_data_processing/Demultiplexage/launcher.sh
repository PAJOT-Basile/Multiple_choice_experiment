#!/bin/bash
#SBATCH -A sexisol
#SBATCH --partition fast
#SBATCH --cpus-per-task 100
#SBATCH --mem 100GB


module load snakemake stacks
snakemake -s ./Demultiplex_all.snk --cores 100 --rerun-incomplete --retries 3
module unload snakemake stacks