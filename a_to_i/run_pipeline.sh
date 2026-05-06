#!/usr/bin/env bash
# ============================================================
# A-to-I Editing Pipeline — master submission script
#
# Usage (from this directory):
#   bash run_pipeline.sh [options]
#
# Options:
#   --steps  1,2,3,4,5,6  comma-separated steps to run (default: 1,2,3,4,5,6)
#            1 = separate celltypes
#            2 = reditools on celltypes
#            3 = reditools on whole BAMs
#            4 = filter AGTC
#            5 = intersect with loci BED
#            6 = QC plots (frequency + coverage per cell type)
#   --bam-dir  DIR        override BAM_DIR from config
#   --dry-run             print sbatch commands without submitting
# ============================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"

source "${HERE}/config.sh"

# ---- Argument parsing ----
STEPS="1,2,3,4,5,6"
DRY_RUN=0
CUSTOM_BAM_DIR=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --steps)   STEPS="$2";         shift 2 ;;
    --bam-dir) CUSTOM_BAM_DIR="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1;          shift ;;
    *) echo "[ERROR] Unknown option: $1"; exit 1 ;;
  esac
done

[[ -n "$CUSTOM_BAM_DIR" ]] && BAM_DIR="$CUSTOM_BAM_DIR"

mkdir -p "$LOG_DIR"

run_sbatch() {
  if [[ "$DRY_RUN" -eq 1 ]]; then
    echo "[DRY-RUN] sbatch $*"
    echo "DRY_RUN_JOBID"
  else
    sbatch "$@" | awk '{print $NF}'
  fi
}

step_active() { echo ",$STEPS," | grep -q ",$1,"; }

# ============================================================
# Discover BAMs and build shared list
# ============================================================
echo "[INFO] Scanning BAM_DIR: ${BAM_DIR}"
[[ -d "$BAM_DIR" ]] || { echo "[FATAL] BAM_DIR not found: $BAM_DIR"; exit 1; }

mapfile -t BAMS < <(find "$BAM_DIR" -maxdepth 1 -name "*.bam" ! -name "*.bai" | sort -V)
N_BAMS=${#BAMS[@]}
(( N_BAMS > 0 )) || { echo "[FATAL] No *.bam files in $BAM_DIR"; exit 1; }

BAM_LIST="${LOG_DIR}/bam_list.txt"
printf '%s\n' "${BAMS[@]}" > "$BAM_LIST"
echo "[INFO] Found $N_BAMS BAM(s) -> $BAM_LIST"

# ============================================================
# Step 1: Separate BAMs by cell type
# ============================================================
JID1=""
if step_active 1; then
  echo "[SUBMIT] Step 1: separate celltypes (array 1-${N_BAMS})"
  JID1=$(run_sbatch \
    --job-name="a2i_separate" \
    --partition="${SLURM_PARTITION}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${SLURM_CPUS_SEPARATE}" \
    --mem="${SLURM_MEM_SEPARATE}" \
    --time="${SLURM_TIME_SEPARATE}" \
    --array="1-${N_BAMS}" \
    --mail-type=END,FAIL \
    --mail-user="${SLURM_MAIL}" \
    --chdir="${SCRIPTS_DIR}" \
    --output="${LOG_DIR}/01_separate_%A_%a.out" \
    --error="${LOG_DIR}/01_separate_%A_%a.err" \
    "${SCRIPTS_DIR}/01_separate_celltypes.sbatch")
  echo "[INFO] Step 1 job: $JID1"
fi

# ============================================================
# Step 2: REDItools on cell-type BAMs
# ============================================================
JID2=""
if step_active 2; then
  DEP2=""
  [[ -n "$JID1" ]] && DEP2="--dependency=afterok:${JID1}"
  echo "[SUBMIT] Step 2: reditools on celltypes (array 1-${N_BAMS})"
  JID2=$(run_sbatch \
    --job-name="a2i_redi_celltype" \
    --partition="${SLURM_PARTITION}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${SLURM_CPUS_REDITOOLS}" \
    --mem="${SLURM_MEM_REDITOOLS}" \
    --time="${SLURM_TIME_REDITOOLS}" \
    --array="1-${N_BAMS}" \
    --mail-type=END,FAIL \
    --mail-user="${SLURM_MAIL}" \
    --chdir="${SCRIPTS_DIR}" \
    --output="${LOG_DIR}/02_redi_celltype_%A_%a.out" \
    --error="${LOG_DIR}/02_redi_celltype_%A_%a.err" \
    ${DEP2} \
    "${SCRIPTS_DIR}/02_reditools_celltype.sbatch")
  echo "[INFO] Step 2 job: $JID2"
fi

# ============================================================
# Step 3: REDItools on whole BAMs (optional, runs in parallel with step 2)
# ============================================================
JID3=""
if step_active 3; then
  echo "[SUBMIT] Step 3: reditools on whole BAMs (array 1-${N_BAMS})"
  JID3=$(run_sbatch \
    --job-name="a2i_redi" \
    --partition="${SLURM_PARTITION}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${SLURM_CPUS_REDITOOLS}" \
    --mem="${SLURM_MEM_REDITOOLS}" \
    --time="${SLURM_TIME_REDITOOLS}" \
    --array="1-${N_BAMS}" \
    --mail-type=END,FAIL \
    --mail-user="${SLURM_MAIL}" \
    --chdir="${SCRIPTS_DIR}" \
    --output="${LOG_DIR}/03_redi_%A_%a.out" \
    --error="${LOG_DIR}/03_redi_%A_%a.err" \
    "${SCRIPTS_DIR}/03_reditools.sbatch")
  echo "[INFO] Step 3 job: $JID3"
fi

# ============================================================
# Step 4: Filter AGTC (depends on step 2 and/or 3)
# ============================================================
JID4=""
if step_active 4; then
  DEP4_PARTS=()
  [[ -n "$JID2" ]] && DEP4_PARTS+=("$JID2")
  [[ -n "$JID3" ]] && DEP4_PARTS+=("$JID3")
  DEP4=""
  if [[ ${#DEP4_PARTS[@]} -gt 0 ]]; then
    DEP4="--dependency=afterok:$(IFS=':'; echo "${DEP4_PARTS[*]}")"
  fi
  # Array size: 1 per BED file; over-provision and let extras exit cleanly
  MAX_ARRAY=200
  echo "[SUBMIT] Step 4: filter AGTC (array 1-${MAX_ARRAY})"
  JID4=$(run_sbatch \
    --job-name="a2i_filter" \
    --partition="${SLURM_PARTITION}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${SLURM_CPUS_FILTER}" \
    --mem="${SLURM_MEM_FILTER}" \
    --time="${SLURM_TIME_FILTER}" \
    --array="1-${MAX_ARRAY}" \
    --mail-type=END,FAIL \
    --mail-user="${SLURM_MAIL}" \
    --chdir="${SCRIPTS_DIR}" \
    --output="${LOG_DIR}/04_filter_%A_%a.out" \
    --error="${LOG_DIR}/04_filter_%A_%a.err" \
    ${DEP4} \
    "${SCRIPTS_DIR}/04_filter_agtc.sbatch")
  echo "[INFO] Step 4 job: $JID4"
fi

# ============================================================
# Step 5: Intersect with loci BED (TE or any regions)
# ============================================================
if step_active 5; then
  if [[ -z "${LOCI_BED:-}" ]]; then
    echo "[SKIP] Step 5: LOCI_BED is empty in config.sh — skipping intersection"
  else
    DEP5=""
    [[ -n "$JID4" ]] && DEP5="--dependency=afterok:${JID4}"
    echo "[SUBMIT] Step 5: intersect with loci BED"
    JID5=$(run_sbatch \
      --job-name="a2i_intersect" \
      --partition="${SLURM_PARTITION}" \
      --nodes=1 --ntasks=1 \
      --cpus-per-task="${SLURM_CPUS_INTERSECT}" \
      --mem="${SLURM_MEM_INTERSECT}" \
      --time="${SLURM_TIME_INTERSECT}" \
      --mail-type=END,FAIL \
      --mail-user="${SLURM_MAIL}" \
      --chdir="${SCRIPTS_DIR}" \
      --output="${LOG_DIR}/05_intersect_%j.out" \
      --error="${LOG_DIR}/05_intersect_%j.err" \
      ${DEP5} \
      "${SCRIPTS_DIR}/05_intersect_te.sbatch")
    echo "[INFO] Step 5 job: $JID5"
  fi
fi

# ============================================================
# Step 6: QC plots — frequency + coverage per cell type
# ============================================================
if step_active 6; then
  DEP6=""
  [[ -n "$JID4" ]] && DEP6="--dependency=afterok:${JID4}"
  echo "[SUBMIT] Step 6: QC plots"
  JID6=$(run_sbatch \
    --job-name="a2i_qc" \
    --partition="${SLURM_PARTITION}" \
    --nodes=1 --ntasks=1 \
    --cpus-per-task="${SLURM_CPUS_INTERSECT}" \
    --mem="${SLURM_MEM_INTERSECT}" \
    --time="4:00:00" \
    --mail-type=END,FAIL \
    --mail-user="${SLURM_MAIL}" \
    --chdir="${SCRIPTS_DIR}" \
    --output="${LOG_DIR}/06_qc_%j.out" \
    --error="${LOG_DIR}/06_qc_%j.err" \
    ${DEP6} \
    "${SCRIPTS_DIR}/06_qc_agtc.sbatch")
  echo "[INFO] Step 6 job: $JID6"
fi

echo ""
echo "========================================"
echo "[OK] Pipeline submitted"
echo "  Steps requested : ${STEPS}"
echo "  BAMs            : ${N_BAMS}"
echo "  BAM list        : ${BAM_LIST}"
echo "  Logs            : ${LOG_DIR}"
echo "  OUT_ROOT        : ${OUT_ROOT}"
echo "========================================"
