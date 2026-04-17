#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=08:00:00
#PBS -N snseq_master6

echo "$(date +%T)"

module purge
module load tools/dev
eval "$(~/miniforge3/bin/conda shell.bash hook)"
source activate /rds/general/user/sg3723/home/anaconda3/envs/dmtcp

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline
MAX_JOBS=46
USER_NAME=${USER:-$(whoami)}
MANIFEST="${WD}/sn_outs/sample_manifest.csv"

cd "${WD}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Missing manifest: ${MANIFEST}" >&2
  exit 1
fi

Rscript -e '
setwd("/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline/sn_outs")
signatures <- character(0)
manifest <- read.csv("sample_manifest.csv", stringsAsFactors = FALSE)
expr_ok <- function(sample) {
  status_path <- file.path("by_samples", sample, "expr_filter_status.txt")
  file.exists(status_path) && trimws(readLines(status_path, warn = FALSE, n = 1)) == "ok"
}
for (sample in manifest$sample) {
  if (!expr_ok(sample)) {
    next
  }
  status_path <- file.path("by_samples", sample, "infercna_status.txt")
  if (!file.exists(status_path) || trimws(readLines(status_path, warn = FALSE, n = 1)) != "ok") {
    next
  }
  sig_path <- file.path("by_samples", sample, paste0(sample, "_signatures.rds"))
  if (!file.exists(sig_path)) {
    next
  }
  sig <- readRDS(sig_path)
  signatures <- unique(c(signatures, sig))
}
write.table(signatures, file = "cancer_signatures.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
'

submitted_samples=()
skipped_samples=()

while read -r sample; do
  [[ -n "${sample}" ]] || continue
  while [[ $(qstat | grep "${USER_NAME}" | wc -l) -gt ${MAX_JOBS} ]]; do
    sleep 180
  done

  expr_status_file="sn_outs/by_samples/${sample}/expr_filter_status.txt"
  expr_status_value=""
  if [[ -f "${expr_status_file}" ]]; then
    expr_status_value=$(head -n 1 "${expr_status_file}" | tr -d '\r')
  fi

  status_file="sn_outs/by_samples/${sample}/infercna_status.txt"
  status_value=""
  if [[ -f "${status_file}" ]]; then
    status_value=$(head -n 1 "${status_file}" | tr -d '\r')
  fi

  if [[ "${expr_status_value}" == "ok" && "${status_value}" == "ok" ]]; then
    qsub -v sample="${sample}" -N "${sample}" 6_Malignancy.sh
    submitted_samples+=("${sample}")
  else
    skipped_samples+=("${sample}:expr=${expr_status_value:-missing_status};infercna=${status_value:-missing_status}")
  fi
done < <(tail -n +2 "${MANIFEST}" | cut -d, -f1 | tr -d '"')

echo
echo "Jobs submitted (infercna_status=ok): ${#submitted_samples[@]}"
((${#submitted_samples[@]})) && printf '  %s\n' "${submitted_samples[@]}"

echo
echo "Skipped samples: ${#skipped_samples[@]}"
((${#skipped_samples[@]})) && printf '  %s\n' "${skipped_samples[@]}"

echo
echo "$(date +%T)"
