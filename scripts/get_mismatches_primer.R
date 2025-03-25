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

primers_expand <- fread(args[3], header=F, sep="\t")

if (ncol(primers_expand) == 5) {
	names(primers_expand) <- c("ori_primer", "primer_set", "fwd_primer", "rev_primer", "probe")
} else {
	names(primers_expand) <- c("ori_primer", "primer_set", "fwd_primer", "rev_primer")
}

# Extract the header of the first sequence
header <- names(aln)[1]  # Get the header of the first sequence

# Split the header by spaces and colon to extract relevant fields
header_fields <- strsplit(header, split='[:| --]+')[[1]]
fwd_ref <- header_fields[6]
rev_ref <- header_fields[7]
fwd_length <- nchar(fwd_ref)
rev_length <- nchar(rev_ref)
primerset <- header_fields[4]

# Extract forward and reverse subsequences from the alignment
fwd <- DNAStringSet(aln, start=1, width=fwd_length)
rev <- reverseComplement(DNAStringSet(aln, start=unique(width(aln))-rev_length+1, width=rev_length))

# Function to compare sequences to reference and create mismatch table
compare_to_reference <- function(sequences, reference, type) {
	nseq <- length(sequences)
	ref_len <- nchar(reference)
	
	# Initialize an empty data frame to store mismatch information
	mismatch_table <- data.frame(primer_set=primerset, 
								primer_seq=reference, 
								position = 1:ref_len, 
								reference = strsplit(reference, "")[[1]],
								stringsAsFactors = FALSE)
	
	tablelist <- list()
	all_sites <- list()
	
	# Loop through each sequence and compare to reference
	for (i in 1:nseq) {
		seq <- as.character(sequences[[i]])  # Convert sequence to character
		mygenome <- strsplit(names(sequences), split=":")[[i]][[1]]
		
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
										ref_length=ref_len,
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
		mutate(mutation = case_when(target!="0" ~ paste0(reference, as.character(position), target), TRUE ~ "NA")) %>% 
		mutate(type=type, mutation = case_when(mutation=="NA" ~ "", TRUE ~ mutation)) %>%
		distinct()

	return(final_result)
}

# Now call the function to compare sequences with a reference (e.g., forward reference)
result_fwd <- compare_to_reference(sequences=fwd, reference=fwd_ref, type="Fwd")
result_rev <- compare_to_reference(sequences=rev, reference=rev_ref, type="Rev")
primer_result <- bind_rows(result_fwd, result_rev)

primer_result <- primer_result %>%
	#mutate(ori_primer=primer_set) %>%
	left_join(select(primers_expand, ori_primer, primer_set), by=c("primer_set")) %>%
	relocate(ori_primer, .after=primer_set) %>%
	relocate(genome, .after=primer_seq)

# Write primer result
if (nrow(primer_result) > 0) {
	fwrite(primer_result, file=output, col.names=T, row.names=F, sep="\t", quote=F)

	# # Convert to gt table
	# gt_table <- gt(primer_result) %>%
	# 	tab_spanner(label = "Primer Mismatches", columns = starts_with("position")) %>%
	# 	tab_spanner(label = "Sequence Mismatches", columns = starts_with("mutation")) 
	
	# # Save gt table as PDF (ensure the webshot library is installed for PDF generation)
	# gt_table %>%
	# 	# The file path for your PDF
	# 	gtsave(filename = paste0(tools::file_path_sans_ext(output), "_result.pdf"), path = ".")
} else {
	file.create(output)
}

