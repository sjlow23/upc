#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(Biostrings)
library(dplyr)
library(seqinr)

args <- commandArgs(trailingOnly=TRUE)

# Convert probe fasta to tabular format
probes <- readDNAStringSet(args[1])
probedf <- data.frame(probe = names(probes),
					ref = as.character(probes))

blast <- fread(args[2], header=F, sep="\t")

process_df <- function(mydf, outfile) {
	names(mydf) <- c("probe", "genome", "positions", "primerset", "slen", "length", 
				  "pident", "nident", "mismatch", "sseq")

	mydf <- mydf %>%
		left_join(probedf, by="probe") %>%
		#group_by(probe, genome, ref) %>% 
		#group_by(genome, primerset) %>%
		group_by(genome, probe, ref) %>%
		filter(pident==max(pident) & nident==max(nident)) %>%
		#break ties
		slice_head(n=1) %>%
		group_by(genome, probe, ref) %>%
		#group_by(genome, primerset, ref) %>%
		summarize(sequence=toString(unique(sseq))) %>%
		distinct() %>%
		ungroup() %>%
		mutate(header = paste0(genome, "--", probe, "--", ref, sep="")) 
		#mutate(header = paste0(genome, "--", primerset, "--", ref, sep=""))
	
	write.fasta(sequences=as.list(mydf$sequence), names=mydf$header, file.out=outfile)
	
}

process_df(blast, outfile=args[3])
