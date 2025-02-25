#!/usr/bin/env Rscript

# Load required libraries
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

# Read collated primer file
primer_result <- fread(args[1], header=T, sep="\t")

output_pergenome <- args[2]
output_permutation_combo <- args[3]
output_perprimer <- args[4]
output_pergenome_mismatches <- args[5]
output_grouped <- args[6]

origenomecount <- as.numeric(args[7])

# Get combinations of primers, genomes, and type
combinations <- expand.grid(
	genome = unique(primer_result$genome),
	ori_primer = unique(primer_result$ori_primer),
	type = unique(primer_result$type)
)

# Number of genomes amplified per primerset (including perfect matches)
amplified_summary <- primer_result %>%
  group_by(ori_primer) %>%
  summarize(count_genomes_amplified=n_distinct(genome),
  			perc_genomes_amplified=round(n_distinct(genome)/origenomecount*100, 2))


# Positions with mismatches to keep in df
keep <- primer_result %>%
	filter(target!=0) %>%
	distinct(ori_primer, position) 

# Limit primers df to only those primers and positions with at least one mismatch
primer_result <- primer_result %>%
	inner_join(keep)

# Per primerset ori and genome combo, number of mismatches, positions of mismatches, 
# and actual mismatches
summary_pergenome <- primer_result %>% 
	mutate(mismatch=case_when(mutation=="" ~ 0, TRUE ~ 1)) %>%
	group_by(genome, ori_primer, position, type) %>% 
	filter(mismatch==min(mismatch)) %>%
	filter(mismatch!=0) %>% 
	ungroup() %>% 
	group_by(genome, ori_primer, primer_seq, type) %>% 
	summarize(count_mismatches=n_distinct(mutation), 
						positions=toString(unique(position)),
						mutations=toString(unique(mutation)),
						binding_sites=toString(unique(binding_site)))

if (nrow(summary_pergenome) > 0) {
	fwrite(summary_pergenome, file=output_pergenome, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_pergenome)
}


# Collapse pergenome output into mutation combinations frequency per primerset
# Alert levels for mismatches: 
# High: mutation at last base
# Medium: mutation at any of last 3 bases, OR at least 2 mutations
# Low: mutation at any other location

summary_permutation_combo <- summary_pergenome %>%
	select(-count_mismatches) %>%
	group_by(ori_primer, primer_seq, binding_sites, type, mutations) %>% 
	summarize(count_genomes_with_mutation = n_distinct(genome), 
	          genomes_with_mutation = toString(unique(genome)),
			  primer_sequence = toString(unique(primer_seq)),
              binding_site = toString(unique(binding_sites))) %>%
	ungroup() %>%
	left_join(amplified_summary, by = "ori_primer") %>%
	mutate(count_db_genomes = origenomecount,
		   perc_with_mutation_amplified = round(count_genomes_with_mutation/count_genomes_amplified*100, 2),
		   perc_with_mutation_total = round(count_genomes_with_mutation/origenomecount*100, 2)) %>%
	mutate(position = str_extract_all(mutations, "(?<=\\D)(\\d+)(?=\\D)")) %>%
	mutate(position = sapply(position, function(x) paste(x, collapse = ", ")),
			count_mismatches = str_count(mutations, pattern=",") + 1) %>%
	mutate(alert = case_when(position == nchar(binding_site) ~ "High",
							 position >= nchar(binding_site) - 2 | count_mismatches >= 2 ~ "Medium",
							TRUE ~ "Low")) %>%
	relocate(genomes_with_mutation, .after = last_col()) %>%
	select(-primer_seq, -binding_sites, -count_mismatches) %>%
  	relocate(count_genomes_with_mutation, .after = binding_site) %>%
	relocate(alert, .after = mutations) %>%
  	arrange(ori_primer, type)


if (nrow(summary_permutation_combo) > 0) {
	fwrite(summary_permutation_combo, file=output_permutation_combo, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_permutation_combo)
}


# Per primerset ori and position combo, mutations
summary_perprimer <- primer_result %>% 
	mutate(mismatch=case_when(mutation == "" ~ 0, TRUE ~ 1)) %>%
	group_by(ori_primer, genome, position, type) %>% 
	filter(mismatch == min(mismatch)) %>%
	filter(mismatch != 0) %>% 
	ungroup() %>%
	group_by(ori_primer, position, mutation, type) %>% 
	summarize(genomes_with_mutation = toString(unique(genome)),
						count_genomes_with_mutation = n_distinct(genome))

if (nrow(summary_perprimer) > 0) {
	fwrite(summary_perprimer, file=output_perprimer, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_perprimer)
}

# Per genome, number of fwd and rev mismatches
# Can't add all possible combinations, as those that fail to amplify are not included in ori primers tsv
# Stats is only for genomes that had amplification
summary_pergenome_mismatches <- primer_result %>% 
	mutate(mismatch = case_when(mutation=="" ~ 0, TRUE ~ 1)) %>%
	group_by(genome, ori_primer, position, type) %>% 
	filter(mismatch == min(mismatch)) %>%
	filter(mismatch != 0) %>% 
	ungroup() %>% 
	group_by(genome, ori_primer, type) %>%
	summarize(mismatches = n_distinct(position))

if (nrow(summary_pergenome_mismatches) > 0) {
	fwrite(summary_pergenome_mismatches, file=output_pergenome_mismatches, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_pergenome_mismatches)
}


# Add in genomes with perfect matches
present <- summary_pergenome_mismatches %>% 
	select(genome, ori_primer, type) %>%
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


















