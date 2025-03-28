#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(knitr)
library(kableExtra)


args <- commandArgs(trailingOnly=TRUE)

primerlist <- fread(args[1], header=F, sep="\t")
probelist <- fread(args[2], header=F, sep="\t")

primer_mismatches <- fread(args[3], header=T, sep="\t")
primer_missing <- fread(args[4], header=T, sep="\t")
probe_mismatches <- fread(args[5], header=T, sep="\t")
probe_missing <- fread(args[6], header=T, sep="\t")

primer_combo_plot <- args[7]
primer_mutation_plot <- args[8]
probe_combo_plot <- args[9]
probe_mutation_plot <- args[10]

primer_pie <- args[11]
probe_pie <- args[12]
status_alluvial <- args[13]
if (file.size(args[14]) > 0) {
	primer_msa_plot <- args[14]
} else {
	primer_msa_plot <- NULL
}

if (file.size(args[15]) > 0) {
	probe_msa_plot <- args[15]
} else {
	probe_msa_plot <- NULL
}

oriprimers <- fread(args[16], header=F, sep="\t")
output_html <- args[17]

names(primerlist) <- c("Primer original", "Primer derived", "Forward primer sequence (5' -> 3')", "Reverse primer sequence (5' -> 3')")
names(probelist) <- c("Primer original", "Probe derived", "Probe sequence (5' -> 3')")
probelist <- probelist %>%
	relocate("Probe derived")



# Find degenerate positions in primers and probes

if (ncol(oriprimers) == 4) {
	names(oriprimers) <- c("primer", "fwd", "rev", "probe")
	probepresent <- "yes"
} else {
	names(oriprimers) <- c("primer", "fwd", "rev")
	probepresent <- "no"
}


find_degenerate_positions <- function(df, column_name) {
  # Define the allowed characters set (A, C, G, T), case-insensitive
  allowed_characters <- c("A", "C", "G", "T")
  
  # Initialize a vector to store the invalid positions for each row
  degenerate_positions_vector <- character(nrow(df))
  
  # Loop through each row in the data frame
  for (i in 1:nrow(df)) {
	# Get the sequence from the specified column
	sequence <- df[[column_name]][i]
	
	# Initialize a vector to hold invalid positions for the current sequence
	degenerate_positions <- integer(0)
	
	# Check each character in the sequence
	for (j in seq_along(strsplit(sequence, NULL)[[1]])) {
	  char <- toupper(substr(sequence, j, j))  # Convert to uppercase and get character
	  # If the character is not in the allowed set, record its position (1-based index)
	  if (!(char %in% allowed_characters)) {
		degenerate_positions <- c(degenerate_positions, j)
	  }
	}
	
	# Convert invalid positions to a comma-separated string and store it in the result vector
	degenerate_positions_vector[i] <- paste(degenerate_positions, collapse = ", ")
  }
  
  return(degenerate_positions_vector)
}


# Join derived primers with degenerate positions
oriprimers$fwd_degenerate <- find_degenerate_positions(oriprimers, "fwd")
oriprimers$rev_degenerate <- find_degenerate_positions(oriprimers, "rev")
if (probepresent == "yes") {
	oriprimers$probe_degenerate <- find_degenerate_positions(oriprimers, "probe")
}


primerlist <- primerlist %>%
	left_join(oriprimers, by = c("Primer original" = "primer")) %>%
	relocate(fwd_degenerate, .after = "Forward primer sequence (5' -> 3')") %>%
	relocate(rev_degenerate, .after = "Reverse primer sequence (5' -> 3')") %>%
	select(`Primer original`, `Primer derived`, `Forward primer sequence (5' -> 3')`, `Reverse primer sequence (5' -> 3')`, fwd_degenerate, rev_degenerate)

probelist <- probelist %>%
	left_join(oriprimers, by = c("Primer original" = "primer")) %>%
	relocate(probe_degenerate, .after = "Probe sequence (5' -> 3')") %>%
	select(`Probe derived`, `Primer original`, `Probe sequence (5' -> 3')`, probe_degenerate)


# Function to highlight degenerate positions
highlight_degenerate <- function(sequence, positions) {
  pos <- as.numeric(strsplit(positions, ",")[[1]])
  chars <- strsplit(sequence, "")[[1]]
  
  # Loop through and highlight text
  for (p in pos) {
	chars[p] <- paste0("<span style='background-color:rgb(172, 255, 127); color: black;'>", chars[p], "</span>")
  }
  return(paste0(chars, collapse = ""))
}

# Original primer table with highlighted degenerate positions
oriprimers$fwd <- mapply(highlight_degenerate, oriprimers$fwd, oriprimers$fwd_degenerate)
oriprimers$rev <- mapply(highlight_degenerate, oriprimers$rev, oriprimers$rev_degenerate)
oriprimers <- oriprimers %>%
	rename(`Forward primer sequence (5' -> 3')` = fwd,
			`Reverse primer sequence (5' -> 3')` = rev)

if (probepresent == "yes") {
	oriprimers$probe <- mapply(highlight_degenerate, oriprimers$probe, oriprimers$probe_degenerate)
	oriprimers <- oriprimers %>%
		rename(`Probe sequence (5' -> 3')` = probe)
}

oriprimers <- oriprimers %>%
	select(-ends_with("degenerate"))

# Highlight degenerate positions in derived primer and probe lists
primerlist$`Forward primer sequence (5' -> 3')` <- mapply(highlight_degenerate, primerlist$`Forward primer sequence (5' -> 3')`, primerlist$fwd_degenerate)
primerlist$`Reverse primer sequence (5' -> 3')` <- mapply(highlight_degenerate, primerlist$`Reverse primer sequence (5' -> 3')`, primerlist$rev_degenerate)
if (probepresent == "yes") {
	probelist$`Probe sequence (5' -> 3')` <- mapply(highlight_degenerate, probelist$`Probe sequence (5' -> 3')`, probelist$probe_degenerate)
}

primerlist <- primerlist %>% 
	rename(`Fwd degenerate position(s)` = fwd_degenerate,
			`Rev degenerate position(s)` = rev_degenerate)
probelist <- probelist %>% 
	rename(`Probe degenerate position(s)` = probe_degenerate)



# Function to highlight characters in the binding_site column
highlight_mismatch <- function(binding_site, positions) {
	pos <- as.numeric(strsplit(positions, ",")[[1]])
	chars <- strsplit(binding_site, "")[[1]]

	# Loop through and highlight text
	for (p in pos) {
		chars[p] <- paste0("<span style='background-color: #FF7f7f; color: black;'>", chars[p], "</span>")
	}
	return(paste0(chars, collapse = ""))
}


########################################################################################################################
# Apply the function to highlight text
primer_mismatches$binding_site <- mapply(highlight_mismatch, primer_mismatches$binding_site, primer_mismatches$position)
probe_mismatches$binding_site <- mapply(highlight_mismatch, probe_mismatches$binding_site, probe_mismatches$position)


# Subset relevant columns if needed
primer_mismatches <- primer_mismatches %>%
		select(ori_primer, primer_set, type, mutations, alert, binding_site, count_db_genomes,
		count_genomes_amplified, perc_genomes_amplified, count_genomes_with_mutation, perc_with_mutation_amplified,
		perc_with_mutation_total, genomes_with_mutation) %>%
	rename(Primer = ori_primer,
		`Primer derived` = primer_set,
		Type = type,
		"Mismatch(es)" = mutations,
		Alert = alert,
		"Primer binding site" = binding_site,
		"Number of genomes in database" = count_db_genomes,
		"Number of genomes amplified" = count_genomes_amplified,
		"Percentage of genomes amplified" = perc_genomes_amplified,
		"Number of genomes with >=1 mismatch(es)" = count_genomes_with_mutation,
		"Percentage of genomes with >= 1 mismatch(es)" = perc_with_mutation_amplified,
		"Percentage of genomes with mismatch(es) (of total in db)" = perc_with_mutation_total,
		"Genome(s) with mismatch(es)" = genomes_with_mutation)

primer_mismatches_fwd <- primer_mismatches %>%
	filter(Type == "Fwd") %>%
	select(-Type) %>%
	arrange(desc(`Percentage of genomes with >= 1 mismatch(es)`))
primer_mismatches_rev <- primer_mismatches %>%
	filter(Type == "Rev") %>%
	select(-Type) %>%
	arrange(desc(`Percentage of genomes with >= 1 mismatch(es)`))

probe_mismatches <- probe_mismatches %>%
	select(ori_probe, mutations, alert, binding_site, count_db_genomes, 
		count_genomes_present, perc_genomes_present, count_genomes_with_mutation, perc_with_mutation_present, 
		perc_with_mutation_total, genomes_with_mutation) %>%
	rename(Probe = ori_probe, 
		"Mismatch(es)" = mutations, 
		Alert = alert, 
		"Probe binding site" = binding_site, 
		"Number of genomes in database" = count_db_genomes, 
		"Number of genomes with probe binding site" = count_genomes_present, 
		"Percentage of genomes with probe binding site" = perc_genomes_present,
		"Number of genomes with >=1 mismatch(es)" = count_genomes_with_mutation, 
		"Percentage of genomes with >=1 mismatch(es)" = perc_with_mutation_present, 
		"Percentage of genomes with mismatch(es) (of total in db)" = perc_with_mutation_total, 
		"Genome(s) with mismatch(es)" = genomes_with_mutation) %>%
	arrange(desc(`Percentage of genomes with >=1 mismatch(es)`))

########################################################################################################################
# Generate missing table
primer_missing <- primer_missing %>%
	rename(Primer = ori_primer,
			Type = type,
			"Number of genomes amplified" = count_genomes_amplified,
			"Number of genomes not amplified" = count_genomes_missing,
			"Genome(s) not amplified" = genomes_missing)

probe_missing <- probe_missing %>%
	rename(Probe = ori_primer,
			"Number of genomes with probe binding site" = count_genomes_present,
			"Number of genomes without probe binding site" = count_genomes_missing,
			"Genome(s) without probe binding site" = genomes_missing)

########################################################################################################################
# Create the table with kable and HTML formatting
title_primer_fwd <- "<h5><strong>Table 3:</strong> Mismatches observed in the forward primer binding region</h5>"
table_mismatch_primer_fwd_html <- kable(primer_mismatches_fwd, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
	column_spec(2, width = "5em") %>%  
  	column_spec(3, width = "8em") %>%  
  	column_spec(4, width = "4em") %>%
	column_spec(5, width = "10em") %>%  
  	column_spec(6, width = "4em") %>%  
  	column_spec(7, width = "4em") %>%
	column_spec(8, width = "4em") %>%  
  	column_spec(9, width = "4em") %>%  
  	column_spec(10, width = "4em") %>%
	column_spec(11, width = "4em") %>%  
  	column_spec(12, width = "30em") 

title_primer_rev <- "<h5><strong>Table 4:</strong> Mismatches observed in the reverse primer binding region</h5>"
table_mismatch_primer_rev_html <- kable(primer_mismatches_rev, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
	column_spec(2, width = "5em") %>%  
  	column_spec(3, width = "8em") %>%  
  	column_spec(4, width = "4em") %>%
	column_spec(5, width = "10em") %>%  
  	column_spec(6, width = "4em") %>%  
  	column_spec(7, width = "4em") %>%
	column_spec(8, width = "4em") %>%  
  	column_spec(9, width = "4em") %>%  
  	column_spec(10, width = "4em") %>%
	column_spec(11, width = "4em") %>%  
  	column_spec(12, width = "30em") 

title_probe <- "<h5><strong>Table 6:</strong> Mismatches observed in the probe binding region</h5>"
table_mismatch_probe_html <- kable(probe_mismatches, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
  	column_spec(2, width = "8em") %>%  
  	column_spec(3, width = "4em") %>%
	column_spec(4, width = "10em") %>%  
  	column_spec(5, width = "4em") %>%  
  	column_spec(6, width = "4em") %>%
	column_spec(7, width = "4em") %>%  
  	column_spec(8, width = "4em") %>%  
  	column_spec(9, width = "4em") %>%
	column_spec(10, width = "4em") %>%  
  	column_spec(11, width = "30em") 

title_primer_missing <- "<h5><strong>Table 7:</strong> Genomes not amplified by primers</h5>"
table_primer_missing <- kable(primer_missing, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
  	column_spec(2, width = "3em") %>%  
  	column_spec(3, width = "5em") %>%
	column_spec(4, width = "5em") %>%  
  	column_spec(5, width = "30em") 

title_probe_missing <- "<h5><strong>Table 8:</strong> Genomes without probe binding sites</h5>"
table_probe_missing <- kable(probe_missing, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
  	column_spec(2, width = "7em") %>%  
	column_spec(3, width = "7em") %>%  
  	column_spec(4, width = "20em")
	

########################################################################################################################
# Embed plots into the HTML
primerpieplot <- paste("<img src='", basename(primer_pie), "' alt='Primer <i>in silico</i> amplification results' width='650'>", 
			   "<h5><strong>Figure 1:</strong> Primer <i>in silico</i> amplification results</h5>", sep = "")

probepieplot <- paste("<img src='", basename(probe_pie), "' alt='Probe <i>in silico</i> amplification results' width='650'>", 
			   "<h5><strong>Figure 2:</strong> Probe <i>in silico</i> amplification results</h5>", sep = "")

primercomboplot <- paste("<img src='", basename(primer_combo_plot), "' alt='Prevalence of mismatch combinations (primer)' width='1600'>", 
			   			"<h5><strong>Figure 3:</strong> Prevalence of mismatch combinations (primer)</h5>", sep = "")

primermutationplot <- paste("<img src='", basename(primer_mutation_plot), "' alt='Proportion of genomes with mismatches in primer (individual positions)' width='1300'>", 
							"<h5><strong>Figure 4:</strong> Proportion of genomes with mismatches in primer (individual positions)</h5>", sep = "")

probecomboplot <- paste("<img src='", basename(probe_combo_plot), "' alt='Prevalence of mismatch combinations (probe)' width='1600'>", 
						"<h5><strong>Figure 5:</strong> Prevalence of mismatch combinations (probe)</h5>", sep = "")

probemutationplot <- paste("<img src='", basename(probe_mutation_plot), "' alt='Proportion of genomes with mismatches in probe (individual positions)' width='1300'>", 
			   "<h5><strong>Figure 6:</strong> Proportion of genomes with mismatches in probe (individual positions)</h5>", sep = "")

alluvialplot <- paste("<img src='", basename(status_alluvial), "' alt='Genome status based on primer and probe binding' width='500'>", 
			   "<h5><strong>Figure 7:</strong> Genome profiles of primer and probe binding sites</h5>", sep = "")

if (!is.null(primer_msa_plot)) {
	primermsaplot <- paste("<img src='", basename(primer_msa_plot), "' alt='MSA of primers' width='1200'>",
				"<h5><strong>Figure 8:</strong> MSA of primers for non-amplified genome(s)</h5>", sep = "")
}

if (!is.null(probe_msa_plot)) {
	probemsaplot <- paste("<img src='", basename(probe_msa_plot), "' alt='MSA of probes' width='600'>",
				"<h5><strong>Figure 9:</strong> MSA of probes not found in genome(s)</h5>", sep = "")
}


########################################################################################################################
# Supplementary
title_oriprimers <- "<h5><strong>Table 1:</strong> Original list of primers</h5>"
table_oriprimers <- kable(oriprimers, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
  	column_spec(2, width = "5em") %>%  
  	column_spec(3, width = "5em") %>%
	column_spec(4, width = "20em")

title_primerlist <- "<h5><strong>Table 2:</strong> List of primers evaluated</h5>"
table_primerlist <- kable(primerlist, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
  	column_spec(2, width = "3em") %>%  
  	column_spec(3, width = "5em") %>%
	column_spec(4, width = "5em") %>%  
  	column_spec(5, width = "5em") %>%
	column_spec(6, width = "30em")

title_probelist <- "<h5><strong>Table 5:</strong> List of probes evaluated</h5>"
table_probelist <- kable(probelist, "html", escape = FALSE) %>%
	kable_styling(bootstrap_options = c("striped", "hover"), html_font = "Monospace") %>%
	column_spec(1, width = "5em") %>%  
  	column_spec(2, width = "5em") %>%  
  	column_spec(3, width = "5em") %>%
	column_spec(4, width = "20em")


# Combine tables with their titles and the plot
combined_html <- paste(title_oriprimers, table_oriprimers, "<br>", "<br>",
					   title_primerlist, table_primerlist, "<br>", "<br>",
					   title_primer_fwd, table_mismatch_primer_fwd_html, "<br>", "<br>",
					   title_primer_rev, table_mismatch_primer_rev_html, "<br>", "<br>",
					   title_probelist, table_probelist, "<br>", "<br>",
					   title_probe, table_mismatch_probe_html, "<br>", "<br>",
					   title_primer_missing, table_primer_missing, "<br>", "<br>",
					   title_probe_missing, table_probe_missing, "<br>", "<br>",
					   primerpieplot, "<br>", "<br>",
					   probepieplot, "<br>", "<br>",
					   primercomboplot, "<br>", "<br>",
					   primermutationplot, "<br>", "<br>",
					   probecomboplot, "<br>", "<br>",
					   probemutationplot, "<br>", "<br>",
					   alluvialplot, "<br>", "<br>",
					   #primermsaplot, "<br>", "<br>",
					   #probemsaplot, "<br>", "<br>",
					   sep = "")

# Allow for optional MSA plot
if (exists("primermsaplot")) {
  combined_html <- paste(combined_html, primermsaplot, "<br>", "<br>", sep = "")
}
if (exists("probemsaplot")) {
  combined_html <- paste(combined_html, probemsaplot, "<br>", "<br>", sep = "")
}


# Save the HTML file
save_kable(combined_html, file = output_html)
