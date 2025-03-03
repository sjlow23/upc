#!/usr/bin/env Rscript

# Load required libraries
library(tidyr)
library(dplyr)
library(data.table)


args <- commandArgs(trailingOnly=TRUE)

aln <- fread(args[2], header=FALSE, sep="\t")
output <- args[3]
names(aln) <- c("representative", "primer_set", "size", "fwd", "rev", "alignment")

if (file.exists(args[1])) {
	# if duplicates present
	lookup <- fread(args[1], header=F, sep="\t")
	names(lookup) <- c("genome", "representative")
	
	lookup <- lookup %>%
	separate_rows(genome, sep=", ") %>%
	full_join(aln, by = "representative") %>%
	mutate(genome = case_when(is.na(genome) ~ representative,
								TRUE ~ genome))
	
	} else {
	cat("No duplicates found\n")
	lookup <- aln %>% 
		mutate(genome = representative) %>%
		relocate(genome, .after=representative)	
	}

fwrite(lookup, output, sep=",", quote=F, row.names=F, col.names=F)


