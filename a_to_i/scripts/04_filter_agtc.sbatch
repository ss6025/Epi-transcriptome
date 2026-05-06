#!/bin/bash
# Step 4 — Filter REDItools output to A-to-I editing sites (AG / TC substitutions)
# Applies MIN_COV and MIN_FREQ thresholds from config.sh
# Processes all *_RNA_redi.bed files from both celltype and whole-BAM reditools output.
# One file per array task (task 1 builds the file list; extras exit cleanly).

#SBATCH --job-name=a2i_filter
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=75G
#SBATCH --time=12:00:00
#SBATCH --array=1-200

set -euo pipefail
export LC_ALL=C

SCRIPTS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
source "${SCRIPTS_DIR}/../config.sh"

mkdir -p "$LOG_DIR" "$FILTERED_AGTC_OUT"

RUNLOG="${LOG_DIR}/04_filter_runtime_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}.log"
exec > >(tee -a "$RUNLOG") 2>&1

echo "============================================================"
echo "[BOOT] $(date)"
echo "[JOB ] ${SLURM_JOB_ID:-NA}[${SLURM_ARRAY_TASK_ID:-NA}]"
echo "============================================================"

LIST="${FILTERED_AGTC_OUT}/inputs.list"
TID="${SLURM_ARRAY_TASK_ID:-1}"

# Task 1 builds the list of all input BED files to process
if [[ "$TID" -eq 1 ]]; then
  {
    # Cell-type reditools output (individual cell types, not the combined file)
    find "${REDITOOLS_CELLTYPE_OUT}" -maxdepth 1 -type f \
      -name "*_RNA_redi.bed" ! -name "*_all_celltypes_*" 2>/dev/null || true
    # Whole-BAM reditools output
    find "${REDITOOLS_OUT}" -maxdepth 1 -type f \
      -name "*_RNA_redi.bed" 2>/dev/null || true
  } | sort -V > "${LIST}.tmp"
  mv "${LIST}.tmp" "$LIST"
  echo "[INFO] Found $(wc -l < "$LIST") input BED files"
else
  # Wait for task 1 to write the list
  for _ in $(seq 1 60); do
    [[ -s "$LIST" ]] && break
    sleep 2
  done
fi

[[ -s "$LIST" ]] || { echo "[FATAL] Input list not ready: $LIST"; exit 2; }

N=$(wc -l < "$LIST")
(( TID <= N )) || { echo "[INFO] Task $TID > N=$N; nothing to do."; exit 0; }

IN=$(sed -n "${TID}p" "$LIST")
bn=$(basename "$IN")
SAMPLE="${bn%_RNA_redi.bed}"

OUT_BED="${FILTERED_AGTC_OUT}/${SAMPLE}.AGTC.bed"

echo "[INFO] IN       = $IN"
echo "[INFO] OUT      = $OUT_BED"
echo "[INFO] MIN_COV  = ${MIN_COV}"
echo "[INFO] MIN_FREQ = ${MIN_FREQ}"

# ---- AWK filter: keep AG/TC sites passing coverage and frequency thresholds ----
awk -v SAMPLE="$SAMPLE" -v MINCOV="${MIN_COV}" -v MINFREQ="${MIN_FREQ}" -F'\t' '
BEGIN {
  OFS="\t"
  chr_i=0; pos_i=0; strand_i=0; cov_i=0; freq_i=0; subs_i=0; celltype_i=0
}

function lc(s) { gsub(/[[:space:]]+/, "", s); return tolower(s) }

NR==1 {
  for (i=1; i<=NF; i++) {
    key = lc($i)
    if (key=="region"   || key=="chr"        || key=="chrom" || key=="chromosome") chr_i=i
    if (key=="position" || key=="pos")                                              pos_i=i
    if (key=="strand")                                                              strand_i=i
    if (key=="coverage-q30" || key=="coverageq30" || key=="coverage")             cov_i=i
    if (key=="frequency" || key=="freq")                                           freq_i=i
    if (key=="allsubs"   || key=="subs"      || key=="sub")                        subs_i=i
    if (key=="celltype"  || key=="cell_type" || key=="cell")                       celltype_i=i
  }
  # Fall back to positional layout if header detection failed
  if (cov_i==0 || freq_i==0 || subs_i==0 || pos_i==0 || chr_i==0) {
    chr_i=1; pos_i=2; strand_i=4; cov_i=5; subs_i=8; freq_i=9; celltype_i=NF
  }
  next
}

{
  chr     = $(chr_i)
  pos     = $(pos_i)
  if (pos !~ /^[0-9]+$/) next

  allsubs = $(subs_i)
  if (allsubs != "AG" && allsubs != "TC") next

  cov  = (cov_i  > 0 ? $(cov_i)+0  : 0)
  freq = (freq_i > 0 ? $(freq_i)+0 : 0)

  if (cov  < MINCOV)  next
  if (freq < MINFREQ) next

  start    = pos - 1
  end      = pos
  strand   = (strand_i   > 0 ? $(strand_i)   : ".")
  celltype = (celltype_i > 0 ? $(celltype_i) : ".")
  name     = SAMPLE "|" chr ":" pos "|" allsubs

  print chr, start, end, name, 0, strand, allsubs, cov, freq, celltype
}
' "$IN" \
| { echo -e "Chrom\tStart\tEnd\tName\tScore\tStrand\tAllSubs\tCoverage\tFrequency\tCellType"; cat; } \
> "$OUT_BED"

echo "[OK] Wrote $(wc -l < "$OUT_BED") rows -> $OUT_BED"
echo "[DONE] $(date)"
