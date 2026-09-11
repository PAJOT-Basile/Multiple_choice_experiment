#!/bin/bash
#SBATCH -A sexisol
#SBATCH --partition fast
#SBATCH --cpus-per-task 1
#SBATCH --mem 64GB

infolder="/shared/projects/sexisol/archive/25_ddRAD_novogene_demultiplexed"

for sample_id in $(ls ${infolder}/*_L[1-9].[12].fq.gz | cut -d"/" -f7 | cut -d"_" -f1 | sort | uniq); do
   initial_nb_lines=$(zcat ${infolder}/${sample_id}_L[0-9].[12].fq.gz | wc -l)
   second_nb_lines=$(zcat ${infolder}/${sample_id}.[12].fq.gz | wc -l)
   if [[ "${initial_nb_lines}" == "${second_nb_lines}" ]]; then
      echo -e "${sample_id}\tOk"
   else
      echo -e "${sample_id}\tAAAAAAAAAh!"
   fi
done