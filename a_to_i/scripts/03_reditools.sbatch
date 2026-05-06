#!/bin/bash
# Step 3 — Run REDItools on whole BAMs (one per array task)
# Outputs: {sample}_RNA_redi.bed under REDITOOLS_OUT
# All paths come from config.sh

#SBATCH --job-name=a2i_redi
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=72:00:00

set -euo pipefail
umask 002
export LC_ALL=C
export PYTHONNOUSERSITE=1

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPTS_DIR}/../config.sh"

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

mkdir -p "$LOG_DIR"
RUNLOG="${LOG_DIR}/03_redi_runtime_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}.log"
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

# ---- Get BAM for this task ----
BAM_LIST="${LOG_DIR}/bam_list.txt"
[[ -s "$BAM_LIST" ]] || { echo "[FATAL] BAM list not found: $BAM_LIST"; exit 1; }

N=$(wc -l < "$BAM_LIST")
TID="${SLURM_ARRAY_TASK_ID:-1}"
(( TID >= 1 && TID <= N )) || { echo "[INFO] Task $TID > N=$N; nothing to do."; exit 0; }

BAM=$(sed -n "${TID}p" "$BAM_LIST")
[[ -s "$BAM" ]] || { echo "[FATAL] BAM not found: $BAM"; exit 1; }

b="$(basename "$BAM")"
SAMPLE="${b%.bam}"

mkdir -p "${REDITOOLS_OUT}"
WORK_DIR="${REDITOOLS_OUT}/_work/${SAMPLE}"
mkdir -p "$WORK_DIR"

WHOLE_BED="${REDITOOLS_OUT}/${SAMPLE}_RNA_redi.bed"

echo "[INFO] BAM    = $BAM"
echo "[INFO] SAMPLE = $SAMPLE"
echo "[INFO] OUT    = $WHOLE_BED"

# ---- Build REF-matching BAM (filter contigs to those in reference) ----
REF_FAI="${SEQ_REF}.fai"
[[ -s "$REF_FAI" ]] || "$SAM" faidx "$SEQ_REF"

REF_CONTIGS="${WORK_DIR}/ref.contigs.txt"
BAM_CONTIGS="${WORK_DIR}/bam.contigs.txt"
KEEP_CONTIGS="${WORK_DIR}/${SAMPLE}.keep.refcontigs.txt"

cut -f1 "$REF_FAI" | sort -u > "$REF_CONTIGS"
"$SAM" idxstats "$BAM" | cut -f1 | grep -v '^\*$' | sort -u > "$BAM_CONTIGS"
comm -12 "$REF_CONTIGS" "$BAM_CONTIGS" > "$KEEP_CONTIGS"

NKEEP=$(wc -l < "$KEEP_CONTIGS" | awk '{print $1}')
echo "[INFO] REF∩BAM contigs = $NKEEP"
(( NKEEP >= 5 )) || { echo "[FATAL] Too few shared contigs ($NKEEP). Check genome build."; exit 1; }

REFMATCH_BAM="${WORK_DIR}/${SAMPLE}.refmatch.sorted.bam"
TMP_PFX="${WORK_DIR}/sort.${SLURM_JOB_ID:-0}_${TID}"

"$SAM" view -h "$BAM" \
| awk -v FS="\t" -v OFS="\t" -v keep="$KEEP_CONTIGS" '
  BEGIN{ while((getline<keep)>0){ ok[$1]=1 } }
  /^@/{
    if($1=="@SQ"){
      sn=""
      for(i=1;i<=NF;i++){ if($i~/^SN:/){ sn=substr($i,4); break } }
      if(sn=="" || ok[sn]) print
    } else print
    next
  }
  { r=$3; rn=$7
    if(!ok[r]) next
    if(rn!="=" && rn!="*" && !ok[rn]) next
    print
  }' \
| "$SAM" view -b - \
| "$SAM" sort -@ "$OMP_NUM_THREADS" -T "$TMP_PFX" -o "$REFMATCH_BAM" -

"$SAM" index "$REFMATCH_BAM"

# ---- Run REDItools ----
REDI_STDERR="${WORK_DIR}/${SAMPLE}.reditools.stderr.log"
echo "[INFO] Running REDItools..."
set +e
"$PY" "$REDITOOLS_PATH" \
    -f "$REFMATCH_BAM" \
    -r "$SEQ_REF" \
    -o "$WHOLE_BED" \
    -t "$OMP_NUM_THREADS" \
    2> "$REDI_STDERR"
RC=$?
set -e

echo "[INFO] REDItools exit_code=$RC"
if (( RC != 0 )); then
  echo "[ERROR] REDItools failed — see $REDI_STDERR"
  exit 1
fi
[[ -s "$WHOLE_BED" ]] || echo "[WARN] Empty output BED for $SAMPLE"

echo "[DONE] $(date) $SAMPLE"
