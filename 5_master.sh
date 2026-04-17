#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=08:00:00
#PBS -N snseq_master5

echo "$(date +%T)"

module purge
module load tools/dev

WD=/rds/general/ephemeral/project/tumourheterogeneity1/ephemeral/snSeq_Pipeline
MAX_JOBS=46
USER_NAME=${USER:-$(whoami)}
MANIFEST="${WD}/sn_outs/sample_manifest.csv"

cd "${WD}"

if [[ ! -f "${MANIFEST}" ]]; then
  echo "Missing manifest: ${MANIFEST}" >&2
  exit 1
fi

submitted_samples=()
skipped_samples=()

while read -r sample; do
  [[ -n "${sample}" ]] || continue
  while [[ $(qstat | grep "${USER_NAME}" | wc -l) -gt ${MAX_JOBS} ]]; do
    sleep 180
  done

  status_file="sn_outs/by_samples/${sample}/expr_filter_status.txt"
  status_value=""
  if [[ -f "${status_file}" ]]; then
    status_value=$(head -n 1 "${status_file}" | tr -d '\r')
  fi

  if [[ "${status_value}" == "ok" ]]; then
    qsub -v sample="${sample}" -N "${sample}" 5_InferCNA.sh
    submitted_samples+=("${sample}")
  else
    skipped_samples+=("${sample}:${status_value:-missing_status}")
  fi
done < <(tail -n +2 "${MANIFEST}" | cut -d, -f1 | tr -d '"')

echo
echo "Jobs submitted (expr_filter_status=ok): ${#submitted_samples[@]}"
((${#submitted_samples[@]})) && printf '  %s\n' "${submitted_samples[@]}"

echo
echo "Skipped samples: ${#skipped_samples[@]}"
((${#skipped_samples[@]})) && printf '  %s\n' "${skipped_samples[@]}"

echo
echo "$(date +%T)"
