#!/bin/bash
#SBATCH --job-name=redi_WHOLE_DIR
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=rotembc1@mskcc.org
#SBATCH --chdir=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts
#SBATCH --output=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs/redi_WHOLE_DIR_%A_%a.out
#SBATCH --error=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs/redi_WHOLE_DIR_%A_%a.err

set -eo pipefail
umask 002
export LC_ALL=C
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-16}"

# -------------------- LOGGING --------------------
LOG_DIR="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs"
mkdir -p "$LOG_DIR"
RUNLOG="${LOG_DIR}/redi_WHOLE_DIR_runtime_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"
exec > >(tee -a "$RUNLOG") 2>&1
trap 'echo "[FATAL] line=$LINENO exit=$? cmd=$BASH_COMMAND"; exit 1' ERR

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] job_id=${SLURM_JOB_ID:-NA}"
echo "[TASK] array_task_id=${SLURM_ARRAY_TASK_ID:-NA}"
echo "[NODE] $(hostname)"
echo "[PWD ] $(pwd)"
echo "============================================================"

# -------------------- ENV (conda) --------------------
eval "$(conda shell.bash hook)"
conda activate reditools2

PY="$(command -v python)"
SAM="$(command -v samtools)"
echo "[INFO] python   = $PY"
echo "[INFO] samtools = $SAM"

python - <<'PY'
import pysam
print("[SYSTEM] PYSAM VERSION", pysam.__version__)
print("[SYSTEM] PYSAM PATH", pysam.__path__)
PY

# -------------------- PATHS --------------------
BAM_DIR="/data1/greenbab/users/suns3/MDanderson/bulkRNA_tissue/ALN/sorted_bams"

mapfile -t BAMS < <(ls -1 "${BAM_DIR}"/*MM38.sorted.bam 2>/dev/null | sort -V)

N=${#BAMS[@]}
echo "[INFO] N_BAMS = $N"
printf "[INFO] BAM_LIST:\n%s\n" "${BAMS[@]}"

(( N > 0 )) || { echo "[FATAL] no *.bam found in $BAM_DIR"; exit 1; }

IDX=$((SLURM_ARRAY_TASK_ID - 1))
if (( IDX < 0 || IDX >= N )); then
  echo "[FATAL] task id ${SLURM_ARRAY_TASK_ID} out of range (1..${N})"
  exit 1
fi

BAM="${BAMS[$IDX]}"
b="$(basename "$BAM")"
SAMPLE="${b%.bam}"   # yields R787-1-2762_MM38 etc.

PIPE_BASE="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline"
OUT_BASE="${PIPE_BASE}/tissue"
OUT_WHOLE="${OUT_BASE}/Whole_2"
WORK="${OUT_BASE}/_work_whole_dir"
mkdir -p "$OUT_WHOLE" "$WORK"

SEQ_REF="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/seq_ref_mm10/mm10.fa"
REF="${SEQ_REF}"

RT="/data1/greenbab/users/rotembc1/pipeline2.0/REDItools2/src/cineca/reditools.py"

echo "[INFO] BAM_DIR  = $BAM_DIR"
echo "[INFO] REF      = $REF"
echo "[INFO] REDItools= $RT"
echo "[INFO] OUT_WHOLE = $OUT_WHOLE"
echo "[INFO] WORK      = $WORK"

[[ -d "$BAM_DIR" ]] || { echo "[FATAL] missing BAM_DIR: $BAM_DIR"; exit 1; }
[[ -s "$REF" ]] || { echo "[FATAL] missing REF: $REF"; exit 1; }
[[ -s "${REF}.fai" ]] || "$SAM" faidx "$REF"
[[ -f "$RT" ]] || { echo "[FATAL] missing REDItools script: $RT"; exit 1; }
# -------------------- build REF-matching BAM --------------------
SAMPLE_WORK="${WORK}/${SAMPLE}"
mkdir -p "$SAMPLE_WORK"

REF_FAI="${REF}.fai"
REF_CONTIGS="${SAMPLE_WORK}/ref.contigs.txt"
BAM_CONTIGS="${SAMPLE_WORK}/bam.contigs.txt"
KEEP_CONTIGS="${SAMPLE_WORK}/${SAMPLE}.keep.refcontigs.txt"

cut -f1 "$REF_FAI" | sort -u > "$REF_CONTIGS"
"$SAM" idxstats "$BAM" | cut -f1 | grep -v '^\*$' | sort -u > "$BAM_CONTIGS"
comm -12 "$REF_CONTIGS" "$BAM_CONTIGS" > "$KEEP_CONTIGS"

NKEEP=$(wc -l < "$KEEP_CONTIGS" | awk '{print $1}')
echo "[INFO] REF∩BAM contigs = $NKEEP"
if (( NKEEP < 5 )); then
  echo "[FATAL] Too few contigs overlap between BAM and REF. Check build/species."
  echo "        head REF contigs:"; head -5 "$REF_CONTIGS" || true
  echo "        head BAM contigs:"; head -5 "$BAM_CONTIGS" || true
  exit 1
fi

REFMATCH_SORT="${SAMPLE_WORK}/${SAMPLE}.refmatch.sorted.bam"
# Filter BAM to contigs shared with reference, then sort+index (avoids contig mismatch errors).

echo "[INFO] Building REF-matching BAM -> $REFMATCH_SORT"
"$SAM" view -h "$BAM" \
| awk -v FS="\t" -v OFS="\t" -v keep="$KEEP_CONTIGS" '
BEGIN{ while((getline<keep)>0){ ok[$1]=1 } }
/^@/{
  if($1=="@SQ"){
    sn=""
    for(i=1;i<=NF;i++){ if($i ~ /^SN:/){ sn=substr($i,4); break } }
    if(sn=="" || ok[sn]) print
  } else print
  next
}
{
  r=$3; rn=$7
  if(!ok[r]) next
  if(rn!="=" && rn!="*" && !ok[rn]) next
  print
}' \
| "$SAM" view -b - \
| "$SAM" sort -@ "$OMP_NUM_THREADS" -o "$REFMATCH_SORT" -

"$SAM" index "$REFMATCH_SORT"

# -------------------- run REDItools WHOLE --------------------
WHOLE_BED="${OUT_WHOLE}/${SAMPLE}_RNA_redi.bed"
REDI_STDERR="${SAMPLE_WORK}/${SAMPLE}.reditools.stderr.log"

echo "[INFO] Running REDItools WHOLE -> $WHOLE_BED"
set +e
"$PY" "$RT" -f "$REFMATCH_SORT" -r "$REF" -o "$WHOLE_BED" -t "${OMP_NUM_THREADS}" 2> "$REDI_STDERR"
rc=$?
set -e
echo "[INFO] REDItools exit_code=$rc"

if [[ ! -s "$WHOLE_BED" ]]; then
  echo "[WARN] REDItools produced empty output for $SAMPLE."
  echo "[WARN] See stderr: $REDI_STDERR"
  exit 0
fi

echo "[DONE] $(date) wrote:"
echo "  $WHOLE_BED"
