#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)


args <- commandArgs(trailingOnly=TRUE)

primers <- fread(args[1], header=T, sep="\t")
probes <- fread(args[2], header=T, sep="\t")

primerpie_output <- args[3]
probepie_output <- args[4]

primers <- primers %>%
		group_by(ori_primer, status) %>%
		summarize(value = n_distinct(genome))

# probes <- probes %>%
#     group_by(probe, status) %>%
#     summarize(value = n_distinct(genome))


# This assumes one degenerate probe per primer set
probes_grouped <- probes %>% 
	ungroup() %>%
	group_by(genome, ori_primer) %>%
	summarise(status_all = toString(status)) %>%
	mutate(best_status = case_when(grepl("perfect match", status_all) ~ "Probe perfect match",
																 grepl("with mismatch", status_all) ~ "Probe with mismatch",
																 grepl("Not present", status_all) ~ "Not present",
																 TRUE ~ NA)) %>%
	ungroup() %>%
	group_by(ori_primer, best_status) %>%
	summarize(value = n_distinct(genome)) 

 
# Compute the position of labels
data_primers <- primers %>% 
	arrange(desc(status), .by_group=T) %>%
	group_by(ori_primer) %>%
	mutate(prop = value / sum(value) *100) %>%
	mutate(ypos = cumsum(prop)- 0.5*prop) %>%
	rename(`Genome status` = status)

data_probes <- probes_grouped %>% 
	arrange(desc(best_status), .by_group=T) %>%
	group_by(ori_primer) %>%
	mutate(prop = value / sum(value) *100) %>%
	mutate(ypos = cumsum(prop)- 0.5*prop) %>%
	rename(`Genome status` = best_status)


primer_colors <- c("#999999" = "Not amplified",
				   "#009E73" = "Primer perfect match",
				   "#E69F00" = "Primer with mismatch")
				   
probe_colors <- c("#999999" = "Not present",
				   "#009E73" = "Probe perfect match",
				   "#E69F00" = "Probe with mismatch")
				   

# Basic piechart
primerpie <- ggplot(data_primers, aes(x="", y=prop, fill=`Genome status`)) +
		geom_bar(stat="identity", width=1, color="white") +
		coord_polar("y", start=0) +
		theme_void() + 
		theme(legend.position="bottom") +
		geom_text_repel(aes(y = ypos, label = paste0(round(prop, 1), "%")), color = "black", size=4) +
		scale_fill_manual(values=names(primer_colors), breaks = primer_colors) +
		facet_wrap(.~ori_primer) +
		theme(legend.text = element_text(size=12),
						legend.title = element_text(size=12),
						strip.text = element_text(size=12))

probepie <- ggplot(data_probes, aes(x="", y=prop, fill=`Genome status`)) +
		geom_bar(stat="identity", width=1, color="white") +
		coord_polar("y", start=0) +
		theme_void() + 
		theme(legend.position="bottom") +
		geom_text_repel(aes(y = ypos, label = paste0(round(prop, 1), "%")), color = "black", size=4) +
		scale_fill_manual(values=names(probe_colors), breaks = probe_colors) +
		facet_wrap(.~ori_primer) +
		theme(legend.text = element_text(size=12),
						legend.title = element_text(size=12),
						strip.text = element_text(size=12))


ggsave(primerpie_output, primerpie, width=10, height=6, bg="white", dpi=800)
ggsave(probepie_output, probepie, width=10, height=6, bg="white", dpi=800)