#!/usr/bin/env Rscript

# Load required libraries
library(Biostrings)
library(dplyr)
library(tidyr)
library(data.table)


args <- commandArgs(trailingOnly=TRUE)

# Read alignment file
aln <- readDNAStringSet(args[1])
output <- args[2]
probes_expand <- fread(args[3], header=F, sep="\t")
names(probes_expand) <- c("ori_primer", "probe_set", "probe")

# Extract the header of the first sequence
header <- names(aln)[1]  # Get the header of the first sequence

# Split the header by double dashes extract relevant fields
header_fields <- strsplit(header, split='--')[[1]]
probe_ref <- header_fields[3]
probe_length <- nchar(probe_ref)
probeset <- header_fields[2]

probe <- DNAStringSet(aln, start=1, width=probe_length)

# Function to compare sequences to reference and create mismatch table
compare_to_reference <- function(sequences, reference, type) {
	nseq <- length(sequences)
	ref_len <- nchar(reference)
	
	# Initialize an empty data frame to store mismatch information
	mismatch_table <- data.frame(probe_set=probeset, 
								probe_seq=reference, 
								position = 1:ref_len, 
								reference = strsplit(reference, "")[[1]],
								stringsAsFactors = FALSE)
	
	tablelist <- list()
	all_sites <- list()
	
	# Loop through each sequence and compare to reference
	for (i in 1:nseq) {
		seq <- as.character(sequences[[i]])  # Convert sequence to character
		mygenome <- strsplit(names(sequences), split="--")[[i]][[1]]
		
		# Identify mismatches
		mismatches <- sapply(1:ref_len, function(pos) {
			if (substring(seq, pos, pos) != substring(reference, pos, pos)) {
				return(substring(seq, pos, pos))  # Return the mismatch
			} else {
				return(as.character(0))  # No mismatch
			}
		})
		
		# Extract the binding site for the current sequence (target sequence)
		target_seq <- seq  # This is the sequence being compared
		
		# Create a data frame for binding site information
		bindingsite_data <- data.frame(genome=mygenome, 
									binding_site=target_seq,
									ref_sequence=reference,
									ref_length=probe_length,
									stringsAsFactors = FALSE)
		
		# Store mismatch results and binding sites
		tablelist[[mygenome]] <- mismatches
		all_sites[[i]] <- bindingsite_data
	}
	
	# Combine the mismatch results into a single data frame
	res <- bind_cols(tablelist)
	
	# Bind the binding site data frames together, pivot to join with targetseq
	all_sites_df <- bind_rows(all_sites) 
	
	# Combine the mismatch data with the binding site info
	final_result <- bind_cols(mismatch_table, res) %>%
		pivot_longer(-c(1:4), names_to="genome", values_to="target")

	# Merge the binding sites with the final result
	final_result <- final_result %>%
		left_join(all_sites_df, by="genome") %>%
		#filter(target!="0") %>%
	 mutate(mutation = case_when(target!=0 ~ paste0(reference, position, target), TRUE ~ "NA")) %>% 
	 mutate(type=type, mutation = case_when(mutation=="NA" ~ "", TRUE ~ mutation)) %>%
	 distinct()
	
	
	return(final_result)
}


# Now call the function to compare sequences with a reference (e.g., forward reference)
probe_result <- compare_to_reference(sequences=probe, reference=probe_ref, type="Probe")

probe_result <- probe_result %>%
	left_join(select(probes_expand, ori_primer, probe_set), by="probe_set") %>%
	#mutate(ori_probe=probe_set) %>%
	#separate(probe_set, into=c("ori_probe", "probe"), sep="_", remove=F) %>%
	rename(ori_probe = probe_set) %>%
	relocate(ori_probe) %>%
	relocate(ori_primer, .after = ori_probe) %>%
	relocate(genome, .after = probe_seq)


# Write probe result
if (nrow(probe_result) > 0) {
	fwrite(probe_result, file=output, col.names=T, row.names=F, sep="\t", quote=F)
} else {
	file.create(output)
}