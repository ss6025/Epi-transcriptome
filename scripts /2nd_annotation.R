#!/usr/bin/Rscript
## intall packages 
library(tools)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
intersect_file <- args[1]
full_file <- args[2]
refPath <- args[3]

ref <- read_tsv(refPath)

# Extract only the base filename
output_file_name <- gsub(".txt", "_annotated.txt", basename(full_file))

rm_ref=read_tsv(intersect_file)
colnames(rm_ref) <- c("Region","Position","RepMask_gid","RepMask_strand","IR")
tmp <- left_join(read_tsv(full_file, col_names = F) %>% select(X1,X2,X3,X4,X5,X6,X7,X8,X9) %>% 
    rename(Region = X1, Position = X2, Reference = X3, Strand = X4,`Coverage-q25` =X5, MeanQ = X6, `BaseCount[A,C,G,T]`=X7, 
           AllSubs=X8, Frequency = X9),ref,  
    by = c("Region","Position","Reference")) %>% 
    left_join(rm_ref) %>%
    filter(Region %in% c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10",
                        "chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20",
                        "chr21","chr22","chrX","chrY")) %>%
                        # the SubTypes has been filtered in the previous step only rows start with AG, TC are retained. 
    #filter(AllSubs %in% c("CT","AC","AT","TC","TA","AG","GA","CA","GC","TG","GT","CG")) %>% 
    select(Region,Position,Reference,Strand, IR, REDI_Rep_feat, MeanQ,`Coverage-q25`,`BaseCount[A,C,G,T]`, AllSubs,Frequency, RepMask_gid, RepMask_strand)
# "Strand_ref","REDI_Rep_feat","REDI_REP_id"        "REDI_RefSeq_feat"   "REDI_RefSeq_gid"    "RepMask_gid"        "RepMask_strand" 
tmp[is.na(tmp)] <- "."   # Replace all NAs by .
write_tsv(drop_na(tmp), output_file_name)

