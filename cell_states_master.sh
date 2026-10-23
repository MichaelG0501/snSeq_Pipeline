#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=08:00:00
#PBS -N snseq_cell_states
#PBS -koed
#PBS -o /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/
#PBS -e /rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline/

####################
# Current centred refined noreg route. See analysis/ANALYSIS_MAP.md and
# analysis/methodology/cell_states/cell_states_master_methodology.md.
# Inputs: live sn_outs/cell_states/centred UCell scores and grouping table;
# snSeq_malignant_epi_meta.rds. Outputs: current centred states, adjusted MP
# scores, group maxima, source tables, and cohort figures in the same live tier.
####################

set -euo pipefail

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/project/tumourheterogeneity1/live/snSeq_Pipeline

cd "${WD}"
Rscript analysis/cell_states/assign_states_approach_b_noreg.R

echo "$(date +%T)"
