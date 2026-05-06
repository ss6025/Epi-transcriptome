#!/bin/bash
# Step 6 — QC plots for filtered AGTC sites (frequency + coverage per cell type)
# Runs as a single job over all *.AGTC.bed files in FILTERED_AGTC_OUT.
# All paths come from config.sh.

#SBATCH --job-name=a2i_qc
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=4:00:00

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPTS_DIR}/../config.sh"

mkdir -p "$LOG_DIR"
RUNLOG="${LOG_DIR}/06_qc_runtime_${SLURM_JOB_ID:-0}.log"
exec > >(tee -a "$RUNLOG") 2>&1
trap 'echo "[FATAL] line=$LINENO exit=$? cmd=$BASH_COMMAND"; exit 1' ERR

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] ${SLURM_JOB_ID:-NA}"
echo "[NODE] $(hostname)"
echo "============================================================"

PY="$(command -v python3)" || { echo "[FATAL] python3 not found"; exit 1; }
echo "[INFO] python3 = $PY"

[[ -d "$FILTERED_AGTC_OUT" ]] || { echo "[FATAL] FILTERED_AGTC_OUT not found: $FILTERED_AGTC_OUT"; exit 1; }

N=$(find "$FILTERED_AGTC_OUT" -maxdepth 1 -name "*.AGTC.bed" | wc -l)
echo "[INFO] Found $N AGTC BED file(s) in $FILTERED_AGTC_OUT"
(( N > 0 )) || { echo "[FATAL] No *.AGTC.bed files found"; exit 1; }

QC_OUT="${FILTERED_AGTC_OUT}/qc_agtc"

python3 "${SCRIPTS_DIR}/qc_agtc.py" \
    --input_dir "$FILTERED_AGTC_OUT" \
    --out_dir   "$QC_OUT"

echo "[DONE] $(date) — QC output: $QC_OUT"
