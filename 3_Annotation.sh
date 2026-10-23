#!/bin/bash
#PBS -l select=1:ncpus=2:mem=128gb
#PBS -l walltime=08:00:00
#PBS -N snseq_annot
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/
#PBS -e /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/

set -eo pipefail

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline

cd "${WD}"

Rscript Annotation.R

echo "$(date +%T)"
