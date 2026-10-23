#!/bin/bash
#PBS -l select=1:ncpus=4:mem=128gb
#PBS -l walltime=08:00:00
#PBS -N snseq_gnmf
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/
#PBS -e /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/
set -euo pipefail

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/gnmf

WD=/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline

cd "${WD}"
Rscript analysis/cell_states/prepare_malignant_epithelial_ucell.R

echo "$(date +%T)"
