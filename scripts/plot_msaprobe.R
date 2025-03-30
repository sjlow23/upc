#!/usr/bin/env Rscript

# Load required libraries
library(Biostrings)
library(ggplot2)
#BiocManager::install("ggmsa")
library(ggmsa)

args <- commandArgs(trailingOnly=TRUE)

probealn <- readDNAStringSet(args[1])
probeplot <- args[2]
probe <- args[3]

# Get alignment length
probe_length <- max(nchar(as.character(probealn)))

# Set plot width and height
probe_width <- probe_length * 0.4
probe_height <- ifelse(length(probealn) <=20, length(probealn)*0.8, length(probealn)*0.4)


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


# Get variable positions for probe sequences
probe_highlight <- identify_snps(probealn)

custom_colors_nomismatch <- data.frame(names = c(LETTERS[1:26], "-"), 
						 color = "#FFFFFF", 
						 stringsAsFactors = FALSE)

# Generate plots
plot_msa <- function(myaln, myhighlightpos, myprobe) {
	if (!is.null(myhighlightpos)) {
		myplot <- ggmsa(myaln, 
					char_width=0.5, 
					color="Chemistry_NT", 
					seq_name=T, 
					border="gray",
					posHighligthed = myhighlightpos) + 
		ggtitle(myprobe) +
		theme(axis.text = element_text(size=15),
				plot.title = element_text(size=15, hjust=0.5))
	} else {
		myplot <- ggmsa(myaln, 
					char_width=0.5, 
					custom_color = custom_colors_nomismatch, 
					seq_name=T, 
					border="gray") + 					
		ggtitle(myprobe) +
		theme(axis.text = element_text(size=15),
				plot.title = element_text(size=15, hjust=0.5))
	}
	return(myplot)
}

if (length(probe_highlight) > 0) {
	probe_msaplot <- plot_msa(probealn, myhighlightpos=probe_highlight, myprobe=probe)
} else {
	probe_msaplot <- plot_msa(probealn, myhighlightpos=NULL, myprobe=probe)
}


# Save plots
ggsave(probeplot, probe_msaplot, width=probe_width, height=probe_height, bg="white")
