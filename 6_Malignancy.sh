#!/bin/bash
#PBS -l select=1:ncpus=2:mem=16gb
#PBS -l walltime=01:00:00
#PBS -N snseq_malign

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline

cd "${WD}"

Rscript Malignancy.R "${sample}"

echo "$(date +%T)"
