#!/usr/bin/env bash
# ============================================================
# A-to-I Editing Pipeline — config.sh
# Pipeline: a_to_i
# Source this file from every script: source "${SCRIPTS_DIR}/config.sh"
# ============================================================

PIPELINE_NAME="a_to_i"

# ---- Script directory (absolute path to this pipeline's scripts/) ----
SCRIPTS_DIR="/data1/greenbab/users/rotembc1/pipeline2.0/A_to_I_editing_pipeline/github/scripts"

# ---- Data roots ----
DATA_ROOT="/data1/greenbab/users/rotembc1/tet2_3M_7M_single_cell/GSE264730/tet2_2205_in_vivo/rnaseq_pipeline"
INDEX_ROOT="/data1/greenbab/database/rnaseq-index"
OUT_ROOT="/data1/greenbab/users/rotembc1/tet2_3M_7M_single_cell/GSE264730/tet2_2205_in_vivo/rnaseq_pipeline"
REF_LIST=(MM38)

# Input BAM directory (all *.bam files here will be discovered automatically)
BAM_DIR="${DATA_ROOT}/bams"

# ---- Cell-type metadata ----
# TSV with columns: barcode, mannualAnnot (or manualAnnot)
METADATA="${DATA_ROOT}/combined_markers_noNA.tsv"

# ---- Reference genome ----
# mm10 / GRCm38 — update if using a different build
SEQ_REF="/data1/greenbab/database/rnaseq-index/mm10/mm10.fa"

# ---- REDItools ----
REDITOOLS_PATH="/data1/greenbab/users/rotembc1/pipeline2.0/REDItools2/src/cineca/reditools.py"
REDITOOLS_CONDA_ENV="reditools2"

# ---- Output subdirectories (derived from OUT_ROOT — edit only if needed) ----
CELLTYPE_BAMS="${OUT_ROOT}/celltype_bams"            # per-sample/{celltypes}/
REDITOOLS_CELLTYPE_OUT="${OUT_ROOT}/reditools_celltype_output"
REDITOOLS_OUT="${OUT_ROOT}/reditools_output"          # whole-BAM reditools (not pseudobulk)
FILTERED_AGTC_OUT="${OUT_ROOT}/filtered_agtc"
TE_INTERSECT_OUT="${OUT_ROOT}/te_intersect"
LOG_DIR="${SCRIPTS_DIR}/logs"

# ---- Filtering thresholds (used in 04_filter_agtc) ----
# Set to 0 to disable each filter
MIN_COV=0     # minimum read coverage at a site
MIN_FREQ=0.00    # minimum editing frequency (0–1)

# ---- Loci BED for intersection (step 05) ----
# Can be TE annotations (SQuIRE), known ADAR targets, repeats, etc.
# Set to "" to skip the intersection step entirely.
LOCI_BED="/data1/greenbab/database/te_annotations/mm10_te.bed"

# Whether LOCI_BED coordinates are 1-based (set 1) or 0-based/BED (set 0)
LOCI_ONE_BASED=0

# ---- Sample-level metadata for subsetting (used in 05_intersect_te) ----
# CSV with columns: sample, genotype, condition, treatment  (header required)
# Set to "" to skip subsetting.
SAMPLE_METADATA="${DATA_ROOT}/sample_metadata.csv"

# ---- SLURM defaults ----
SLURM_PARTITION="componc_cpu"
SLURM_MAIL="rotembc1@mskcc.org"

# Resources per step
SLURM_CPUS_SEPARATE=8
SLURM_MEM_SEPARATE="80G"
SLURM_TIME_SEPARATE="24:00:00"

SLURM_CPUS_REDITOOLS=16
SLURM_MEM_REDITOOLS="120G"
SLURM_TIME_REDITOOLS="24:00:00"

SLURM_CPUS_FILTER=1
SLURM_MEM_FILTER="20G"
SLURM_TIME_FILTER="4:00:00"

SLURM_CPUS_INTERSECT=4
SLURM_MEM_INTERSECT="60G"
SLURM_TIME_INTERSECT="12:00:00"
