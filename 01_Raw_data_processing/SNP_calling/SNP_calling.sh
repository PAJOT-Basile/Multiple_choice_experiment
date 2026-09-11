#!/bin/bash

#SBATCH --partition=long
#SBATCH --account=pacobar
#SBATCH --job-name=SNP_calling
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

# Then, we can start executing the script. The first step is to localise the directory in which the scripts
# are localised
HERE="$(pwd)"

# Then, we separate the input configuration file into two independent configuration files that will be used to execute the snakefile
# If it is the first time you run the snakemake or the configuration file has been modified, it will restart this step
RUN="False"
if [ ! -d "${HERE}/Configuration_files" ]; then
    RUN="True"
elif [[ "$(cat ${HERE}/Configuration_files/Date_modif_config.txt)" != "$(stat -c %Y ${HERE}/configuration_file.yaml)" ]]; then
    RUN="True"
fi
# Then, we run the script that allows to separate the configuration file in two to prepare the run for the snakemake
if [[ "${RUN}" = "True" ]]; then
    printf "\rPreparing configuration file: configuration_file.yaml ..."
    "${HERE}/Scripts_snk/Configuration.sh" -w "${HERE}" -c "configuration_file.yaml"
    printf "\rPreparing configuration file: configuration_file.yaml ...              DONE\n"
fi

# Check if the environment files are already created and if not, create them
if [[ ! -d "${HERE}/Configuration_files/envs/" ]]; then
    printf "\rCreating environments ..."
    "${HERE}/Scripts_snk/Create_envs.sh" -f "configuration_file.yaml" -s "snakefile.snk" -c "${HERE}/Configuration_files/envs/"
    "${HERE}/Scripts_snk/Create_envs.sh" -f "configuration_file.yaml" -s "${HERE}/Scripts_snk/Index_ref_genome.snk" -c "${HERE}/Configuration_files/envs/"
    printf "\rCreating environment ...         DONE\n"
fi

# Where to put temporary files
TMPDIR="$(grep "tmp_path" ${HERE}/configuration_file.yaml | cut -f2 -d'"')"
TMP="${TMPDIR}"
TEMP="${TMPDIR}"
mkdir -p "${TMPDIR}"
export TMPDIR TMP TEMP

module load snakemake
# Run the snakemake to index the reference genome
Location_ref_genome=$(grep "Reference_genome" "${HERE}/configuration_file.yaml" | cut -d'"' -f2)
if [ ! -f "${Location_ref_genome}.fai" ] || [ ! -f "${Location_ref_genome}.amb" ]; then
    printf "\rIndexing reference genome"
    snakemake -s "${HERE}/Scripts_snk/Index_ref_genome.snk" --profile "${HERE}/Cluster_profile" --configfile "${HERE}/Configuration_files/Variables_config.yaml" --quiet all
    printf "\rIndexing reference genome ...     DONE\n"
fi

echo "Starting Snakemake execution"
# Run the snakemake file
snakemake -s "${HERE}/snakefile.snk" --profile "${HERE}/Cluster_profile" --configfile "${HERE}/Configuration_files/Variables_config.yaml"


module unload snakemake
echo "Snakemake execution ...            DONE"
