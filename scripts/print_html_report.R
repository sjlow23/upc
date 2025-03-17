#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(knitr)
library(kableExtra)


args <- commandArgs(trailingOnly=TRUE)

primer_mismatches <- fread(args[1], header=T, sep="\t")
primer_missing <- fread(args[2], header=T, sep="\t")
probe_mismatches <- fread(args[3], header=T, sep="\t")
probe_missing <- fread(args[4], header=T, sep="\t")
not_amplified <- fread(args[5], header=T, sep="\t")

output <- args[6]


# Function to highlight characters in the binding_site column
highlight_text <- function(binding_site, positions) {
	pos <- as.numeric(strsplit(positions, ",")[[1]])
	chars <- strsplit(binding_site, "")[[1]]

	# Loop through and highlight text
	for (p in pos) {
		chars[p] <- paste0("<span style='background-color: #FF7f7f; color: black;'>", chars[p], "</span>")
	}
	return(paste0(chars, collapse = ""))
}

# Apply the function to highlight text
mismatches$binding_site <- mapply(highlight_text, mismatches$binding_site, mismatches$position)

# Create the table with kable and HTML formatting
table_html <- kable(mismatches, "html", escape = FALSE) %>%
	#col.names = c("Ori Primer", "Binding Site", "Position")) %>%
	kable_styling(bootstrap_options = c("striped", "hover"))


# Save the HTML file
save_kable(table_html, file = output)