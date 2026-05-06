#!/bin/bash
# Step 5 — Intersect AGTC editing sites with a loci BED (TEs, known sites, etc.)
# LOCI_BED is set in config.sh. Leave it empty to skip this step.
# Runs as a single job (not array) over all AGTC BED files.

#SBATCH --job-name=a2i_intersect
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=60G
#SBATCH --time=12:00:00

set -euo pipefail

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPTS_DIR}/../config.sh"

mkdir -p "$LOG_DIR" "$TE_INTERSECT_OUT"
RUNLOG="${LOG_DIR}/05_intersect_runtime_${SLURM_JOB_ID:-0}.log"
exec > >(tee -a "$RUNLOG") 2>&1

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] ${SLURM_JOB_ID:-NA}"
echo "============================================================"

[[ -n "${LOCI_BED:-}" ]] || { echo "[FATAL] LOCI_BED is empty in config.sh"; exit 1; }
[[ -f "$LOCI_BED"     ]] || { echo "[FATAL] LOCI_BED not found: $LOCI_BED"; exit 1; }

METADATA_ARG=""
[[ -n "${SAMPLE_METADATA:-}" && -f "$SAMPLE_METADATA" ]] && METADATA_ARG="--metadata $SAMPLE_METADATA"

ONE_BASED_ARG=""
[[ "${LOCI_ONE_BASED:-0}" -eq 1 ]] && ONE_BASED_ARG="--loci-one-based"

echo "[INFO] AGTC dir  = $FILTERED_AGTC_OUT"
echo "[INFO] LOCI_BED  = $LOCI_BED"
echo "[INFO] OUT       = $TE_INTERSECT_OUT"

python3 "${SCRIPTS_DIR}/intersect_te.py" \
    --agtc-dir  "$FILTERED_AGTC_OUT" \
    --loci-bed  "$LOCI_BED" \
    --output    "$TE_INTERSECT_OUT" \
    --threads   "${SLURM_CPUS_PER_TASK:-4}" \
    ${ONE_BASED_ARG} \
    ${METADATA_ARG}

echo "[DONE] $(date)"
