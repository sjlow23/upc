#!/usr/bin/env Rscript

# Load required libraries
library(Biostrings)
library(ggplot2)
#BiocManager::install("ggmsa")
library(ggmsa)

args <- commandArgs(trailingOnly=TRUE)

dna_fwd <- readDNAStringSet(args[1])
dna_rev <- readDNAStringSet(args[2])
fwdplot <- args[3]
revplot <- args[4]
primer <- args[5]

# Get alignment length
fwd_length <- max(nchar(as.character(dna_fwd)))
rev_length <- max(nchar(as.character(dna_rev)))

# Set plot width
fwd_width <- fwd_length * 0.4
rev_width <- rev_length * 0.4

# Set plot height
fwd_height <- ifelse(length(dna_fwd) <=20, length(dna_fwd)*0.8, length(dna_fwd)*0.4)
rev_height <- ifelse(length(dna_fwd) <=20, length(dna_rev)*0.8, length(dna_rev)*0.4)

# Identify positions to highlight
identify_snps <- function(dna_set) {
  # Convert aln to character matrix
  dna_matrix <- as.matrix(dna_set)
  
  # Identify positions with SNP
  snp_positions <- apply(dna_matrix, 2, function(col) {
	uniq_values <- unique(col)
	uniq_values <- uniq_values[uniq_values != "-"]
	length(uniq_values) > 1  
  })
  snp_indices <- as.numeric(which(snp_positions))
}

# Get variable positions for forward and reverse DNA sequences
fwd_highlight <- identify_snps(dna_fwd)
rev_highlight <- identify_snps(dna_rev)


# Generate plots
plot_msa <- function(myaln, myhighlightpos, myprimer, mystrand) {
	myplot <- ggmsa(myaln, 
					char_width=0.5, 
					color="Chemistry_NT", 
					seq_name=T, 
					border="gray",
					posHighligthed = myhighlightpos) + 
		ggtitle(paste0(myprimer, ": ", mystrand)) +
		theme(axis.text = element_text(size=15),
				plot.title = element_text(size=15, hjust=0.5))
	return(myplot)
}

fwd_msaplot <- plot_msa(dna_fwd, myhighlightpos=fwd_highlight, myprimer=primer, mystrand="Fwd")
rev_msaplot <- plot_msa(dna_rev, myhighlightpos=rev_highlight, myprimer=primer, mystrand="Rev")

# Save plots
ggsave(fwdplot, fwd_msaplot, width=fwd_width, height=fwd_height, bg="white")
ggsave(revplot, rev_msaplot, width=rev_width, height=rev_height, bg="white")