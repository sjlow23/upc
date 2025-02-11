#!/usr/bin/env Rscript

# Load required libraries
library(data.table)
library(Biostrings)
library(dplyr)
library(seqinr)

args <- commandArgs(trailingOnly=TRUE)


blast_target <- fread(args[1], header=F, sep="\t")
blast_offtarget <- fread(args[2], header=F, sep="\t")


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


process_df(blast_target, outfile=args[3])
process_df(blast_offtarget, outfile=args[4])


#split fasta by probe for target and offtarget
#align, then process similarly to primers

