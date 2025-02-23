#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(Biostrings)
library(dplyr)
library(seqinr)

args <- commandArgs(trailingOnly=TRUE)


blast <- fread(args[1], header=F, sep="\t")


process_df <- function(mydf, outfile) {
	names(mydf) <- c("probe", "genome", "positions", "primerset", "slen", "length", 
				  "pident", "nident", "mismatch", "qseq", "sseq")
	mydf <- mydf %>%
		group_by(probe, genome, sseq) %>% 
		summarize(sequence=toString(unique(qseq))) %>%
		distinct() %>%
		ungroup() %>%
		mutate(header = paste0(genome, "--", probe, "--", sseq, sep="")) 
	
	write.fasta(sequences=as.list(mydf$sequence), names=mydf$header, file.out=outfile)
	
}

process_df(blast, outfile=args[2])
