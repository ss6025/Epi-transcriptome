# Repeat-Aware RNA Editing Quantification Pipeline

A repeat-aware computational pipeline for quantifying RNA editing activity within repetitive elements, enabling systematic analysis of ADAR-mediated A-to-I editing and its impact on innate immune regulation and tumor evolution.

---

## Background & Motivation

RNA editing, particularly adenosine-to-inosine (A-to-I) editing mediated by ADAR enzymes, is a critical post-transcriptional mechanism that modulates double-stranded RNA (dsRNA) sensing and innate immune activation. Repetitive elements (REs)—including Alu, LINE, and ERV families—are the dominant substrates for RNA editing due to their propensity to form dsRNA structures.

However, conventional RNA editing pipelines are known loci dependent or poorly optimized for repetitive regions and frequently exclude multimapping reads or collapse editing signals across loci, obscuring lineage-, context-, and disease-specific editing patterns. This pipeline addresses these limitations by enabling **repeat-aware, quantitative profiling of RNA editing activity**, with a focus on RE-derived dsRNA substrates and their immunological consequences.

---

## Key Features

- Repeat-aware quantification of A-to-I RNA editing across RE families and subfamilies
- Support for Alu-, LINE-, and ERV-centric editing analyses, in addition to gene based analysis
- Identification of editing-enriched inverted repeat (IR) regions
- Compatible with bulk RNA-seq, single-cell or spatial analyses
- Modular design for extension to new editing metrics

---

## Requirements
Install nextflow, java/jdk-11.0.11, singularity/3.7.1
Executor: LSF is default. change profile to local if you need to run it locally.

### Software
- Python ≥ 3.9  
- R ≥ 4.2  
- samtools ≥ 1.15  
- bedtools ≥ 2.30  

### Python Packages
- pandas  
- numpy  
- pysam  
- pybedtools  

### R Packages (optional downstream analysis)
- DESeq2  
- tidyverse  
- ggplot2  

---

## Installation

Clone the repository:

```bash
git clone https://github.com/ss6025/repeat-editing-pipeline.git
cd repeat-editing-pipeline
```

## Inputs
- input_dir: the input directory that contains the bed files (marking the RNA editing positions) from REDITOOLS
- output_dir: directory to store the output (default: results folder under current working directory)
- email_address: email address to be notified when the nextflow run finishes.

## Outputs
outputs will will be saved under the output_dir you specified in the prompt.
- filtered_AI/ : filtered RNA editing txt files
- filtered_AI/annotated/ : annotated files with the corresponding repeat class and families
- filtered_AI/intersect_RM/ : intersected regions with reference
- PDF/ folder: donut.pdf donut_alus.pdf RNA_editing_rate_histogram.pdf
- Global editing tables: GlobalEditing_Freq.csv GlobalEditing_Freq_IRvsNonIR.csv GlobalEditing_Freq_bySubFamilies.csv

Notes: output directory would store the symlinks to original files stored in nextflow working directory. so if you need to copy the actual folder, you can use

```bash
cp -rL <output_dir> <your_destination_folder>
```

## Usage
for Juno users
```bash
module load rna-editing-analysis
rna_editing_analysis.sh 
```
it will prompt you for input_dir, output_dir (default: outputs folder in current working directory), and email address (default is your <msk_user_id>@mskcc.org) then it will kick off nextflow pipeline in the backend. the run.log will store standard out/err. the

for general users
1. clone this repo
2. 
```bash
./rna_editing_analysis.sh 
```
3. answer the prompts and start the run.

## Example
```bash
(base) suns3@conifold:/work/greenbaum/users/lih7/test/random$ rna_editing_analysis.sh
/work/greenbaum/software/rna-editing-analysis
Config file path [/work/greenbaum/software/rna-editing-analysis/nextflow.config]: 
Input directory : /work/greenbaum/users/suns3/MSKCC/iacobuzio/organoids_PDAC/ALN/REDITOOLS_OUTPUT/SINE
Output directory [results]: 
Notify this Email Address [lih7@mskcc.org]: 
Running RNAediting analysis pipeline:
----------------
Input directory: /work/greenbaum/users/suns3/MSKCC/iacobuzio/organoids_PDAC/ALN/REDITOOLS_OUTPUT/SINE
Output directory: /juno/work/greenbaum/users/lih7/test/random/results
----------------

(base) lih7@conifold:/work/greenbaum/users/lih7/test/random$ tree ./
|-- results
|   |-- GlobalEditing_Freq.csv -> /work/greenbaum/software/rna-editing-analysis/work/09/38dd5bd6179531923407f195b4db42/GlobalEditing_Freq.csv
|   |-- GlobalEditing_Freq_IRvsNonIR.csv -> /work/greenbaum/software/rna-editing-analysis/work/09/38dd5bd6179531923407f195b4db42/GlobalEditing_Freq_IRvsNonIR.csv
|   |-- GlobalEditing_Freq_bySubFamilies.csv -> /work/greenbaum/software/rna-editing-analysis/work/09/38dd5bd6179531923407f195b4db42/GlobalEditing_Freq_bySubFamilies.csv
|   |-- PDF
|   |   |-- RNA_editing_rate_histogram.pdf -> /work/greenbaum/software/rna-editing-analysis/work/09/38dd5bd6179531923407f195b4db42/PDF/RNA_editing_rate_histogram.pdf
|   |   |-- donut.pdf -> /work/greenbaum/software/rna-editing-analysis/work/09/38dd5bd6179531923407f195b4db42/PDF/donut.pdf
|   |   `-- donut_alus.pdf -> /work/greenbaum/software/rna-editing-analysis/work/09/38dd5bd6179531923407f195b4db42/PDF/donut_alus.pdf
|   `-- filtered_AI
|       |-- Proj_07323_ABEHJK_s_HT1CSHL_RNA_redi.txt -> /work/greenbaum/software/rna-editing-analysis/work/41/51341e3f7940ce2ff5684d60827061/Proj_07323_ABEHJK_s_HT1CSHL_RNA_redi.txt
|       |-- ...
|   `-- intersect_RM
|       |-- Proj_07323_ABEHJK_s_HT1CSHL_RNA_redi_intersect.txt
|       |-- ...
|   `-- annotated
|       |-- Proj_07323_ABEHJK_s_HT1CSHL_RNA_redi_annotated.txt
|       |-- ...
```

## Contributing
Contributions to this utility are welcome! If you encounter any issues or have suggestions for improvements, please open an issue or submit a pull request on GitHub.

## Authors
Siyu Sun (suns3@mskcc.org)
Mike Li (lih7@mskcc.org)
Camille Rotemberg (rotembc1@mskcc.org)
 
