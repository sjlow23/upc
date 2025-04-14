#!/usr/bin/env Rscript

# Load required libraries
library(tidyr)
library(dplyr)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)

aln <- fread(args[2], header = FALSE, sep = "\t")
output <- args[3]

if (!is.null(args[4])) {
	assembly_lookup <- fread(args[4], header = FALSE, sep= "\t")
}


names(aln) <- c("representative", "primer_set", "size", "fwd", "rev", "alignment")
names(assembly_lookup) <- c("accession", "assembly_accession")

if (file.exists(args[1])) {
	# if duplicates present
	lookup <- fread(args[1], header = F, sep = "\t")
	names(lookup) <- c("genome", "representative")
	
	lookup <- lookup %>%
	separate_rows(genome, sep = ", ") %>%
	full_join(aln, by = "representative") %>%
	mutate(genome = case_when(is.na(genome) ~ representative,
								TRUE ~ genome)) %>%
	select(-representative)
	
} else {
	cat("No duplicates found\n")
	lookup <- aln %>% 
		mutate(genome = representative) %>%
		select(-representative)
		#relocate(genome, .after=representative)	
}

# Convert to assembly accession
if (exists("assembly_lookup")) {
	lookup <- lookup %>% 
		separate(genome, into = c("genome", "position"), sep = ":", remove = T) %>%
		left_join(assembly_lookup, by = c("genome" = "accession")) %>%
		select(-genome) %>%
		mutate(new = paste0(assembly_accession, ":", position)) %>%
		rename(genome = new) %>%
		relocate(genome)
	lookup
}

fwrite(lookup, output, sep=",", quote=F, row.names=F, col.names=F)


