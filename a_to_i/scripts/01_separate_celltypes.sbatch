#!/bin/bash
# Step 1 — Separate BAM by cell type using CB:Z barcodes
# Submitted by run_pipeline.sh with --array=1-N_BAMS
# All paths come from config.sh

#SBATCH --job-name=a2i_separate
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=80G
#SBATCH --time=72:00:00

set -euo pipefail
umask 002
export LC_ALL=C
export PYTHONNOUSERSITE=1

# ---- Load config ----
SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPTS_DIR}/../config.sh"

mkdir -p "$LOG_DIR"
RUNLOG="${LOG_DIR}/01_separate_runtime_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}.log"
exec > >(tee -a "$RUNLOG") 2>&1
trap 'echo "[FATAL] line=$LINENO exit=$? cmd=$BASH_COMMAND"; exit 1' ERR

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]"
echo "[NODE] $(hostname)"
echo "============================================================"

# ---- Validate tools ----
SAM="$(command -v samtools)" || { echo "[FATAL] samtools not found"; exit 1; }
PY="$(command -v python3)"   || { echo "[FATAL] python3 not found";  exit 1; }
echo "[INFO] samtools = $SAM"
echo "[INFO] python3  = $PY"

# ---- Get BAM for this task ----
BAM_LIST="${LOG_DIR}/bam_list.txt"
[[ -s "$BAM_LIST" ]] || { echo "[FATAL] BAM list not found: $BAM_LIST"; exit 1; }

N=$(wc -l < "$BAM_LIST")
TID="${SLURM_ARRAY_TASK_ID:-1}"
(( TID >= 1 && TID <= N )) || { echo "[INFO] Task $TID > N=$N; nothing to do."; exit 0; }

BAM=$(sed -n "${TID}p" "$BAM_LIST")
[[ -s "$BAM" ]] || { echo "[FATAL] BAM not found: $BAM"; exit 1; }

# ---- Validate metadata ----
[[ -f "$METADATA" ]] || { echo "[FATAL] METADATA not found: $METADATA"; exit 1; }

# ---- Derive sample name ----
b="$(basename "$BAM")"
SAMPLE="${b%.bam}"
SAMPLE="${SAMPLE#possorted_genome_}"

OUTPUT_DIR="${CELLTYPE_BAMS}/${SAMPLE}/celltypes"

echo "[INFO] BAM        = $BAM"
echo "[INFO] SAMPLE     = $SAMPLE"
echo "[INFO] METADATA   = $METADATA"
echo "[INFO] OUTPUT_DIR = $OUTPUT_DIR"

# ---- Run separation ----
python3 "${SCRIPTS_DIR}/separate_celltypes.py" \
    --bam     "$BAM" \
    --markers "$METADATA" \
    --output_dir "$OUTPUT_DIR"

echo "[DONE] $(date) $SAMPLE"
