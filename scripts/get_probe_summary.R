#!/usr/bin/env Rscript

# Load required libraries
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

# Read collated primer file
probe_result <- fread(args[1], header=T, sep="\t")

output_pergenome <- args[2]
output_permutation_combo <- args[3]
output_perprobe <- args[4]
output_pergenome_mismatches <- args[5]
output_grouped <- args[6]

origenomecount <- as.numeric(args[7])

# Get combinations of probes, genomes, and type
combinations <- expand.grid(
	genome = unique(probe_result$genome),
	ori_probe = unique(probe_result$ori_probe),
	type = unique(probe_result$type)
)

# Number of genomes amplified per probeset (including perfect matches)
amplified_summary <- probe_result %>%
  group_by(ori_probe) %>%
  summarize(count_genomes_present=n_distinct(genome))


# Positions with mismatches to keep in df
keep <- probe_result %>%
	filter(target!=0) %>%
	distinct(ori_probe, position) 

# Limit probes df to only those probes and positions with at least one mismatch
probe_result <- probe_result %>%
	inner_join(keep)

# Per probeset ori and genome combo, number of mismatches, positions of mismatches, 
# and actual mismatches
summary_pergenome <- probe_result %>% 
	mutate(mismatch=case_when(mutation=="" ~ 0, TRUE ~ 1)) %>%
	group_by(genome, ori_probe, position, type) %>% 
	filter(mismatch==min(mismatch)) %>%
	filter(mismatch!=0) %>% 
	ungroup() %>% 
	group_by(genome, ori_probe, type) %>% 
	summarize(count_mismatches=n_distinct(mutation), 
						positions=toString(unique(position)),
						mutations=toString(unique(mutation)))

if (nrow(summary_pergenome) > 0) {
	fwrite(summary_pergenome, file=output_pergenome, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_pergenome)
}


# Collapse pergenome output into mutation combinations frequency per probeset
summary_permutation_combo <- summary_pergenome %>%
	group_by(ori_probe, type, mutations) %>% 
	summarize(count_genomes_with_mutation=n_distinct(genome), genomes_with_mutation=toString(unique(genome))) %>%
	ungroup() %>%
	left_join(amplified_summary, by="ori_probe") %>%
	mutate(perc_with_mutation_present=round(count_genomes_with_mutation/count_genomes_present*100, 2),
		   perc_with_mutation_total=round(count_genomes_with_mutation/origenomecount*100, 2)) %>%
	mutate(position = str_extract_all(mutations, "(?<=\\D)(\\d+)(?=\\D)")) %>%
	mutate(position = sapply(position, function(x) paste(x, collapse = ", "))) %>%
	relocate(genomes_with_mutation, .after=last_col()) %>%
	arrange(ori_probe, type)
		   

if (nrow(summary_permutation_combo) > 0) {
	fwrite(summary_permutation_combo, file=output_permutation_combo, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_permutation_combo)
}


# Per probeset ori and position combo, mutations
summary_perprobe <- probe_result %>% 
	mutate(mismatch=case_when(mutation=="" ~ 0, TRUE ~ 1)) %>%
	group_by(ori_probe, genome, position, type) %>% 
	filter(mismatch==min(mismatch)) %>%
	filter(mismatch!=0) %>% 
	ungroup() %>%
	group_by(ori_probe, position, mutation, type) %>% 
	summarize(genomes_with_mutation=toString(unique(genome)),
						count_genomes_with_mutation=n_distinct(genome))

if (nrow(summary_perprobe) > 0) {
	fwrite(summary_perprobe, file=output_perprobe, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_perprobe)
}

# Per genome, number of fwd and rev mismatches
# Can't add all possible combinations, as those that fail to amplify are not included in ori primers tsv
# Stats is only for genomes that had amplification
summary_pergenome_mismatches <- probe_result %>% 
	mutate(mismatch=case_when(mutation=="" ~ 0, TRUE ~ 1)) %>%
	group_by(genome, ori_probe, position, type) %>% 
	filter(mismatch==min(mismatch)) %>%
	filter(mismatch!=0) %>% 
	ungroup() %>% 
	group_by(genome, ori_probe, type) %>%
	summarize(mismatches=n_distinct(position))

if (nrow(summary_pergenome_mismatches) > 0) {
	fwrite(summary_pergenome_mismatches, file=output_pergenome_mismatches, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_pergenome_mismatches)
}


# Add in genomes with perfect matches
present <- summary_pergenome_mismatches %>% 
	select(genome, ori_probe, type) %>%
	distinct()
missing_combinations <- combinations %>% 
	anti_join(present)

summary_pergenome_mismatches <- bind_rows(summary_pergenome_mismatches, missing_combinations) %>%
	replace(is.na(.), 0) %>%
	pivot_wider(names_from="type", values_from="mismatches")

if (nrow(summary_pergenome_mismatches) > 0) {
	fwrite(summary_pergenome_mismatches, file=output_grouped, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_grouped)
}



##################################################################################################################















