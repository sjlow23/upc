#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(Biostrings)
library(dplyr)
library(seqinr)
library(tidyr)


args <- commandArgs(trailingOnly=TRUE)

location <- fread(args[1], header=T, sep="\t")
outfile <- args[2]

location <- location %>%
	separate(seqID, into = c("genome", "pos"), sep=":", remove=T) %>%
	rename(probe = patternName, 
			ref = pattern, 
			sequence = matched) %>%
			select(-pos, -strand, -start, -end) %>%
			mutate(header = paste0(genome, "--", probe, "--", ref, sep="")) 

write.fasta(sequences=as.list(location$sequence), names=location$header, file.out=outfile)
