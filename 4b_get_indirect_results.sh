#!/usr/bin/bash -l
#$ -cwd
#$ -N get_indirectresults
#$ -o out.out
#$ -e err.err
#$ -M run7@cdc.gov
#$ -m ae
#$ -t 1-10
#$ -l h_rt=08:00:00

echo "**** New job ****"
date
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

module load gcc
module load R/4.0.4

Rscript 4b_get_indirect_results.R 1000 10 ${SGE_TASK_ID}

module list
echo "End date:"
date
