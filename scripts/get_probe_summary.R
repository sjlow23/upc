#!/usr/bin/env Rscript

# Load required libraries
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)

# Read collated probe file
probe_result <- fread(args[1], header=T, sep="\t")

output_pergenome <- args[2]
output_permutation_combo <- args[3]
output_perprobe <- args[4]
output_pergenome_mismatches <- args[5]
output_grouped <- args[6]
output_missing <- args[7]
output_status <- args[8]
genomelist <- gsub(".fna", "", readLines(args[9]))
origenomecount <- as.numeric(args[10])
probes_expand <- fread(args[11], header=F, sep="\t")
names(probes_expand) <- c("ori_primer", "ori_probe", "probe")
category <- args[12]


# Get combinations of probes, genomes, and type
combinations <- expand.grid(
	genome = unique(probe_result$genome),
	ori_primer = unique(probe_result$ori_primer),
	type = unique(probe_result$type)
)


# Number of genomes amplified per probeset (including perfect matches)
amplified_summary <- probe_result %>%
  group_by(ori_probe) %>%
  summarize(count_genomes_present = n_distinct(genome),
  			perc_genomes_present = round(n_distinct(genome)/origenomecount*100, 2))


# Positions with mismatches to keep in df
keep <- probe_result %>%
	filter(target != 0) %>%
	distinct(ori_probe, position) 

# Limit probes df to only those probes and positions with at least one mismatch
probe_result.subset <- probe_result %>%
	inner_join(keep)

# Per probeset ori and genome combo, number of mismatches, positions of mismatches, 
# and actual mismatches
if (nrow(probe_result.subset) >=1) {
	summary_pergenome <- probe_result.subset %>% 
		mutate(mismatch=case_when(mutation == "" ~ 0, TRUE ~ 1)) %>%
		group_by(genome, ori_probe, position, type) %>% 
		filter(mismatch == min(mismatch)) %>%
		filter(mismatch != 0) %>% 
		ungroup() %>% 
		group_by(genome, ori_probe, probe_seq, type) %>% 
		summarize(count_mismatches = n_distinct(mutation), 
							positions = toString(unique(position)),
							mutations = toString(unique(mutation)),
							binding_sites = toString(unique(binding_site)))
	if (nrow(summary_pergenome) > 0) {
	fwrite(summary_pergenome, file=output_pergenome, col.names=T, row.names=F, sep="\t", quote=F)
	}
} else {
	file.create(output_pergenome)
}


# Get missing genomes
if (category == "target") {
	summary_missing <- probe_result %>%
	group_by(ori_primer) %>%
	summarize(count_genomes_present = n_distinct(genome),
			  count_genomes_missing = origenomecount - count_genomes_present,
			  genomes_missing = toString(setdiff(genomelist, genome)))
} else {
	summary_missing <- probe_result %>%
	group_by(ori_primer) %>%
	summarize(count_genomes_present = n_distinct(genome),
			  count_genomes_missing = origenomecount - count_genomes_present,
			  genomes_present = toString(unique(genome)))
}

if (nrow(summary_missing) > 0) {
	fwrite(summary_missing, file=output_missing, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output_missing)
}



# Collapse pergenome output into mutation combinations frequency per probeset
# Alert levels for probe mismatches: 
# Medium: at least 2 mutations
# Low: mutation at any other location

if (exists("summary_pergenome")) {
	summary_permutation_combo <- summary_pergenome %>%
	left_join(select(probes_expand, ori_probe, ori_primer), by="ori_probe") %>%
	select(-count_mismatches) %>%
	group_by(ori_probe, ori_primer, probe_seq, binding_sites, type, mutations) %>% 
	summarize(count_genomes_with_mutation = n_distinct(genome), 
				genomes_with_mutation = toString(unique(genome)),
				probe_sequence = toString(unique(probe_seq)),
				binding_site = toString(unique(binding_sites))) %>%
	ungroup() %>%
	left_join(amplified_summary, by="ori_probe") %>%
	mutate(count_db_genomes = origenomecount,
		   perc_with_mutation_present = round(count_genomes_with_mutation/count_genomes_present*100, 2),
		   perc_with_mutation_total = round(count_genomes_with_mutation/origenomecount*100, 2)) %>%
	mutate(position = str_extract_all(mutations, "(?<=\\D)(\\d+)(?=\\D)")) %>%
	mutate(position = sapply(position, function(x) paste(x, collapse = ", ")),
		   count_mismatches = str_count(mutations, pattern=",") + 1) %>%
	mutate(alert = case_when(count_mismatches >= 2 ~ "Medium", TRUE ~ "Low")) %>%
	select(-probe_seq, -binding_sites, -count_mismatches) %>%
	relocate(genomes_with_mutation, .after = last_col()) %>%
	relocate(count_genomes_with_mutation, .after = binding_site) %>%
	relocate(alert, .after = mutations) %>%
	arrange(ori_primer, type)

	if (nrow(summary_permutation_combo) > 0) {
	fwrite(summary_permutation_combo, file=output_permutation_combo, col.names=T, row.names=F, sep="\t", quote=F)
	}
} else {
	file.create(output_permutation_combo)
}


# Per probeset ori and position combo, mutations
if (nrow(probe_result.subset) >=1) {
	summary_perprobe <- probe_result.subset %>% 
		mutate(mismatch = case_when(mutation=="" ~ 0, TRUE ~ 1)) %>%
		group_by(ori_probe, genome, position, type) %>% 
		filter(mismatch == min(mismatch)) %>%
		filter(mismatch != 0) %>% 
		ungroup() %>%
		group_by(ori_probe, position, mutation, type) %>% 
		summarize(genomes_with_mutation = toString(unique(genome)),
							count_genomes_with_mutation = n_distinct(genome)) %>%
		left_join(select(probes_expand, ori_probe, ori_primer), by="ori_probe") %>%
		relocate(ori_primer, .after = ori_probe)

	if (nrow(summary_perprobe) > 0) {
		fwrite(summary_perprobe, file=output_perprobe, col.names=T, row.names=F, sep="\t", quote=F)
	}
} else {
	file.create(output_perprobe)
}


# Per genome, number of fwd and rev mismatches
# Can't add all possible combinations, as those that fail to amplify are not included in ori primers tsv
# Stats is only for genomes that had amplification

if (nrow(probe_result.subset) >=1) {
	summary_pergenome_mismatches <- probe_result.subset %>% 
		mutate(mismatch=case_when(mutation == "" ~ 0, TRUE ~ 1)) %>%
		group_by(genome, ori_probe, position, type) %>% 
		filter(mismatch == min(mismatch)) %>%
		filter(mismatch!=0) %>% 
		ungroup() %>% 
		group_by(genome, ori_probe, type) %>%
		summarize(mismatches = n_distinct(position)) %>%
		left_join(select(probes_expand, ori_primer, ori_probe), by="ori_probe")

	if (nrow(summary_pergenome_mismatches) > 0) {
		fwrite(summary_pergenome_mismatches, file=output_pergenome_mismatches, col.names=T, row.names=F, sep="\t", quote=F)
	}
} else {
	file.create(output_pergenome_mismatches)
}


# Add in genomes with perfect matches
if (nrow(probe_result.subset) >=1) {
	present <- summary_pergenome_mismatches %>% 
		select(genome, ori_primer, ori_probe, type) %>%
		distinct()
	missing_combinations <- combinations %>% 
		anti_join(present) %>%
		mutate(ori_probe = as.character(NA)) 
	# Best ori probe match for each genome
	perfectgenomes <- missing_combinations %>% 
		filter(is.na(ori_probe)) %>%
		pull(genome)
	probe_bestmatch <- probe_result %>%
		filter(genome %in% perfectgenomes) %>%
		select(genome, ori_probe) %>%
		distinct() %>%
		dplyr::rename(ori_probe_final = ori_probe)
	summary_pergenome_mismatches_grouped <- bind_rows(summary_pergenome_mismatches, missing_combinations) %>%
		mutate(mismatches = replace(mismatches, is.na(mismatches), 0)) %>%
		left_join(probe_bestmatch, by="genome") %>%
		mutate(ori_probe = case_when(is.na(ori_probe) ~ ori_probe_final, TRUE ~ ori_probe)) %>%
		select(-ori_probe_final)
	fwrite(summary_pergenome_mismatches_grouped, file=output_grouped, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	# All perfect matches
	missing_combinations <- combinations %>% 
		mutate(ori_probe = as.character(NA))
	perfectgenomes <- missing_combinations %>% 
		filter(is.na(ori_probe)) %>%
		pull(genome)
	probe_bestmatch <- probe_result %>%
		filter(genome %in% perfectgenomes) %>%
		select(genome, ori_probe) %>%
		distinct() %>%
		dplyr::rename(ori_probe_final = ori_probe)
	summary_pergenome_mismatches_grouped <- missing_combinations %>%
		mutate(mismatches = 0) %>%
		left_join(probe_bestmatch, by="genome") %>%
		mutate(ori_probe = case_when(is.na(ori_probe) ~ ori_probe_final, TRUE ~ ori_probe)) %>%
		select(-ori_probe_final)
	fwrite(summary_pergenome_mismatches_grouped, file=output_grouped, col.names=T, row.names=F, sep="\t", quote=F)

}



# Summarise per genome status
if (nrow(probe_result.subset) >=1) {
	pergenome_status <- summary_pergenome_mismatches_grouped %>%
		#left_join(select(probes_expand, ori_primer, ori_probe), by="ori_probe") %>%
		mutate(sum_mismatches = sum(mismatches)) %>%
		mutate(status = case_when(sum_mismatches == 0 ~ "Probe perfect match",
								sum_mismatches != 0 ~ "Probe with mismatch")) %>%
		select(genome, ori_primer, status)

	#Unique probes
	uniq_primers <- unique(pergenome_status$ori_primer)

	# List to store missing genomes for each primer set
	missinglist <- list()

	for (primer in uniq_primers) {
		present <- pergenome_status %>%
			filter(ori_primer == primer) %>%
			pull(genome)
		missing <- setdiff(genomelist, present)
		if (length(missing) != 0) {
			missinglist[[primer]] <- data.frame(genome = missing, ori_primer = primer, status = "Not present")
		}
	}

	missingdf <- bind_rows(missinglist)
	pergenome_status <- bind_rows(pergenome_status, missingdf) %>%
		relocate(genome)
	fwrite(pergenome_status, file=output_status, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	# All genomes perfect match
	uniq_primers <- unique(combinations$ori_primer)
	# List to store missing genomes for each primer set
	missinglist <- list()

	for (primer in uniq_primers) {
		missing <- setdiff(genomelist, unique(combinations$genome))
		if (length(missing) != 0) {
			missinglist[[primer]] <- data.frame(genome = missing, ori_probe = "", ori_primer = primer, status = "Not present")
		}
	}

	missingdf <- bind_rows(missinglist)
	
	pergenome_status <- combinations %>%
		mutate(status = "Probe perfect match") %>%
		mutate(ori_probe = "") %>%
		relocate(ori_probe, .after = genome) %>%
		select(-type)
	pergenome_status <- bind_rows(pergenome_status, missingdf)
	fwrite(pergenome_status, file=output_status, col.names=T, row.names=F, sep="\t", quote=F)
}
















