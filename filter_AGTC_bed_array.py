#!/bin/bash
#SBATCH --job-name=bed_AGTC_filt
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24G
#SBATCH --time=2:00:00
#SBATCH --array=1-5
#SBATCH --output=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs/bed_AGTC_filt_%A_%a.out
#SBATCH --error=/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/scripts/logs/bed_AGTC_filt_%A_%a.err
set -euo pipefail

BASE="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline"

IN_DIR="${BASE}/scRNAseq_3M_tissue/Whole"
OUT_DIR="${IN_DIR}/true_bed_full_filtered_AGTC"
mkdir -p "${OUT_DIR}" "${BASE}/scripts/logs"

MIN_COV=25
MIN_FREQ=0.01

LIST="${OUT_DIR}/inputs.list"

if [[ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]]; then
  find "${IN_DIR}" -maxdepth 1 -type f -name "*.bed" -print | sort > "${LIST}.tmp"
  mv "${LIST}.tmp" "${LIST}"
  echo "[INFO] Found $(wc -l < "${LIST}") input .bed files"
else
  for i in $(seq 1 120); do
    [[ -s "${LIST}" ]] && break
    sleep 1
  done
fi

if [[ ! -s "${LIST}" ]]; then
  echo "[FATAL] No input list: ${LIST}"
  exit 2
fi

N=$(wc -l < "${LIST}")
if [[ "${SLURM_ARRAY_TASK_ID}" -gt "${N}" ]]; then
  echo "[INFO] Task ${SLURM_ARRAY_TASK_ID} > N=${N}; nothing to do."
  exit 0
fi

IN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${LIST}")
bn=$(basename "$IN")
sample="${bn%.bed}"

OUT_BED="${OUT_DIR}/${sample}.AGTC.cov${MIN_COV}.freq${MIN_FREQ}.bed"

echo "[INFO] IN      = ${IN}"
echo "[INFO] OUT_BED = ${OUT_BED}"

awk -v MINCOV="${MIN_COV}" -v MINFREQ="${MIN_FREQ}" -v SAMPLE="${sample}" -F'\t' '
BEGIN{
  OFS="\t";
  cov_i=0; freq_i=0; subs_i=0; strand_i=0; pos_i=0; chr_i=0;
  header_done=0;
}
function lc(s){ gsub(/[[:space:]]+/, "", s); return tolower(s) }

NR==1{
  # try to detect header by matching known column names
  for (i=1;i<=NF;i++){
    key=lc($i)
    if (key=="region" || key=="chr" || key=="chrom" || key=="chromosome") chr_i=i
    if (key=="position" || key=="pos") pos_i=i
    if (key=="strand") strand_i=i
    if (key=="coverage-q30" || key=="coverageq30" || key=="coverage") cov_i=i
    if (key=="frequency" || key=="freq") freq_i=i
    if (key=="allsubs" || key=="subs" || key=="sub") subs_i=i
  }

  if (cov_i==0 || freq_i==0 || subs_i==0 || pos_i==0 || chr_i==0){
    # fallback: common REDItools layout
    chr_i=1
    pos_i=2
    strand_i=4
    cov_i=5
    subs_i=8
    freq_i=9
    header_done=0
  } else {
    header_done=1
  }

  if (header_done==1){
    # header exists, skip it
    next
  }
  # else: headerless file; do NOT next — process this line as data below
}

{
  chr = $(chr_i)
  pos = $(pos_i)

  if (pos !~ /^[0-9]+$/) next

  cov  = $(cov_i) + 0
  freq = $(freq_i) + 0
  allsubs = $(subs_i)

  if (cov < MINCOV) next
  if (freq < MINFREQ) next
  if (allsubs != "AG" && allsubs != "TC") next

  start = pos - 1
  end   = pos

  strand = (strand_i>0 ? $(strand_i) : ".")
  name   = SAMPLE "|" chr ":" pos "|" allsubs
  score  = 0

  # BED6 + extras (AllSubs, cov, freq)
  print chr, start, end, name, score, strand, allsubs, cov, freq
}
' "${IN}" > "${OUT_BED}"

echo "[OK] Wrote $(wc -l < "${OUT_BED}") rows -> ${OUT_BED}"
