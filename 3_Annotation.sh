#!/bin/bash
#PBS -l select=1:ncpus=2:mem=128gb
#PBS -l walltime=08:00:00
#PBS -N snseq_annot

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline

cd "${WD}"

Rscript Annotation.R

echo "$(date +%T)"
