#!/bin/bash
#PBS -l select=1:ncpus=4:mem=64gb
#PBS -l walltime=08:00:00
#PBS -N snseq_cell_states
#PBS -koed

####################
# Script: cell_states_master.sh
# Methodology: analysis/methodology/cell_states/cell_states_master_methodology.md
# Map: analysis/ANALYSIS_MAP.md
#
# Description:
#   Run the current snRNA-seq malignant epithelial state workflow:
#   Approach B, noreg, followed by fixed scRef-retained 3CA unresolved relabeling
#   and terminal presentation figures.
#
# Input:
#   sn_outs/snSeq_malignant_epi.rds
#   sn_outs/Metaprogrammes_Results/UCell_nMP19_filtered.rds
#   sn_outs/UCell_3CA_MPs.rds
#   sn_outs/snSeq_malignant_epi_cc_score.rds
#
# Output:
#   sn_outs/Auto_topmp_v2_noreg_states_B.rds
#   sn_outs/Auto_final_states.rds
#   sn_outs/task3_sample_abundance/*
#   sn_outs/task4_unresolved_states/*
#   sn_outs/task6_hybrid_pairwise/*
####################

set -euo pipefail

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline

cd "${WD}"
Rscript analysis/cell_states/assign_states_approach_b_noreg.R
Rscript analysis/cell_states/relabel_unresolved_retained_3ca.R
Rscript analysis/cell_states/plot_state_overall_proportions.R
Rscript analysis/cell_states/plot_state_sample_abundance.R
Rscript analysis/cell_states/plot_state_hybrid_pairwise.R

echo "$(date +%T)"
