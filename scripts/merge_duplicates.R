#!/usr/bin/env Rscript

# Load required libraries
library(tidyr)
library(dplyr)
library(data.table)


args <- commandArgs(trailingOnly = TRUE)

aln <- fread(args[2], header = FALSE, sep = "\t")
output <- args[3]

if (args[4] != "NULL") {
	assembly_lookup <- fread(args[4], header = FALSE, sep= "\t")
    names(assembly_lookup) <- c("accession", "assembly_accession")
}

names(aln) <- c("representative", "primer_set", "size", "fwd", "rev", "alignment")


if (args[1] != "NULL") {
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


