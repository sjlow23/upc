#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(tidyr)
library(data.table)


args <- commandArgs(trailingOnly=TRUE)
# Read alignment file
query <- fread(args[1], header=T, sep="\t")
output <- args[2]
mode <- args[3]

if (mode == "primer") {
    bestmatch <- query %>% 
        ungroup() %>% 
        group_by(ori_primer, primer_set, type, genome) %>% 
        summarize(mismatches=n_distinct(position[target!=0])) %>% 
        ungroup() %>% 
        group_by(genome, ori_primer, primer_set) %>% 
        summarize(sum_mismatch=sum(mismatches)) %>% 
        ungroup() %>% 
        group_by(genome, ori_primer) %>% 
        filter(sum_mismatch==min(sum_mismatch)) %>% 
        slice_head(n=1)

    # Get best combinations to keep
    best_combo <- bestmatch %>% 
        select(genome, ori_primer, primer_set)

    # Inner join with original df
    res <- query %>% 
        inner_join(best_combo)
}

if (mode == "probe") {
    bestmatch <- query %>% 
        ungroup() %>% 
        group_by(ori_probe, ori_primer, type, genome) %>% 
        summarize(mismatches=n_distinct(position[target!=0])) %>% 
        ungroup() %>% 
        group_by(genome, ori_primer, ori_probe) %>% 
        summarize(sum_mismatch=sum(mismatches)) %>% 
        ungroup() %>% 
        group_by(genome, ori_primer) %>% 
        filter(sum_mismatch==min(sum_mismatch)) %>% 
        slice_head(n=1)

    # Get best combinations to keep
    best_combo <- bestmatch %>% 
        select(genome, ori_primer, ori_probe)

    # Inner join with original df
    res <- query %>% 
        inner_join(best_combo)
}


# Write result
if (nrow(res) > 0) {
	fwrite(res, file=output, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output)
}

