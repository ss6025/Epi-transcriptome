#!/bin/bash
# Step 2 — Run REDItools on each cell-type BAM
# Processes one sample's celltypes/ directory per array task
# All paths come from config.sh

#SBATCH --job-name=a2i_redi_celltype
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --time=72:00:00

set -euo pipefail
umask 002
export LC_ALL=C
export PYTHONNOUSERSITE=1

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPTS_DIR}/../config.sh"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$LOG_DIR"
RUNLOG="${LOG_DIR}/02_redi_celltype_runtime_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}.log"
exec > >(tee -a "$RUNLOG") 2>&1
trap 'echo "[FATAL] line=$LINENO exit=$? cmd=$BASH_COMMAND"; exit 1' ERR

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]"
echo "[NODE] $(hostname)"
echo "============================================================"

eval "$(conda shell.bash hook)"
conda activate "${REDITOOLS_CONDA_ENV}"

PY="$(command -v python)"   || PY="$(command -v python3)"
SAM="$(command -v samtools)" || { echo "[FATAL] samtools not found"; exit 1; }
echo "[INFO] python   = $PY"
echo "[INFO] samtools = $SAM"

[[ -s "$SEQ_REF" ]]         || { echo "[FATAL] SEQ_REF not found: $SEQ_REF"; exit 1; }
[[ -f "$REDITOOLS_PATH" ]]  || { echo "[FATAL] REDItools not found: $REDITOOLS_PATH"; exit 1; }

# ---- Discover celltypes directories from step 1 output ----
mapfile -t CELLTYPES_DIRS < <(
  find "$CELLTYPE_BAMS" -mindepth 2 -maxdepth 2 -name "celltypes" -type d | sort -V
)
N=${#CELLTYPES_DIRS[@]}

TID="${SLURM_ARRAY_TASK_ID:-1}"
(( N > 0 ))                              || { echo "[FATAL] No celltypes/ dirs in $CELLTYPE_BAMS"; exit 1; }
(( TID >= 1 && TID <= N ))              || { echo "[INFO] Task $TID > N=$N; nothing to do."; exit 0; }

CELLTYPES_DIR="${CELLTYPES_DIRS[$((TID - 1))]}"

echo "[INFO] celltypes dir = $CELLTYPES_DIR"
echo "[INFO] output dir    = $REDITOOLS_CELLTYPE_OUT"
echo "[INFO] REF           = $SEQ_REF"
echo "[INFO] threads       = $OMP_NUM_THREADS"

mkdir -p "$REDITOOLS_CELLTYPE_OUT"

python3 "${SCRIPTS_DIR}/reditools_on_celltypes.py" \
    --celltypes_dir  "$CELLTYPES_DIR" \
    --output_dir     "$REDITOOLS_CELLTYPE_OUT" \
    --ref            "$SEQ_REF" \
    --reditools_path "$REDITOOLS_PATH" \
    --threads        "$OMP_NUM_THREADS"

echo "[DONE] $(date) $(basename "$(dirname "$CELLTYPES_DIR")")"
