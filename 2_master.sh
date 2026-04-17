#!/bin/bash
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -l walltime=08:00:00
#PBS -N snseq_master2

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

tail -n +2 "${MANIFEST}" | cut -d, -f1 | tr -d '"' | while read -r sample; do
  [[ -n "${sample}" ]] || continue
  while [[ $(qstat | grep "${USER_NAME}" | wc -l) -gt ${MAX_JOBS} ]]; do
    sleep 180
  done
  qsub -v sample="${sample}" -N "${sample}" 2_Clustering.sh
done

echo "$(date +%T)"
