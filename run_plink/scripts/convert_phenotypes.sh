#!/bin/bash

#SBATCH --job-name=makephen
#SBATCH --nodes=1 --mem=64G --time=0-12:00:00

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}


cd ${HOME}/repo/ukbb-instruments/scripts

Rscript convert_phenotypes.r
