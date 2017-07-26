#!/bin/bash

#SBATCH --job-name=runplink
#SBATCH --array=1-200
#SBATCH --nodes=1 --tasks-per-node=28 --time=0-12:00:00

echo "Running on ${HOSTNAME}"
module add R/3.2.3-foss-2016a

if [ -n "${1}" ]; then
  echo "${1}"
  SLURM_ARRAY_TASK_ID=${1}
fi

i=${SLURM_ARRAY_TASK_ID}


cd ${HOME}/repo/ukbb-instruments/scripts

Rscript run_assoc.r ${i} 100
