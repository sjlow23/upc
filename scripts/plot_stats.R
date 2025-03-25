#!/usr/bin/env Rscript

# Load required libraries
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggalluvial)
library(ggmsa)
library(cowplot)
library(wesanderson)


args <- commandArgs(trailingOnly=TRUE)

# Read collated primer file
primer_result <- fread(args[1], header=T, sep="\t")
primer_mutations <- fread(args[2], header=T, sep="\t")
primer_status <- fread(args[3], header=T, sep="\t")
probe_result <- fread(args[4], header=T, sep="\t")
probe_mutations <- fread(args[5], header=T, sep="\t")
probe_grouped <- fread(args[6], header=T, sep="\t")
probe_status <- fread(args[7], header=T, sep="\t")

primeroutput <- args[8]
probeoutput <- args[9]
primer_combo_plot <- args[10]
primer_mutation_plot <- args[11]
probe_combo_plot <- args[12]
probe_mutation_plot <- args[13]
alluvial_plot <- args[14]

mycolors <- c("#FFD700", "#FF0000","#2BE4E9", "#F19CBB", "#ADFF2F", "#2F4F4F", "#BA55D3", "#9C7C38",  
			  "#D8BFD8", "#C0C0C0", "#FF7540", "#107090",  "#D5006A",
			  "#B1A2CA","#EFBE7D","#61F4DE","#FFCACA", "#676767","#FCF5C7","#6B8E23","#87CEFA", 
			  "#E9EC6B",  "#FF4500","#FF1493", "#4D4DFF", "#FF00FF", 
			  "#008000","#008080","#00FF00","#00FF7F", "#acff80", "#d3d3d3", "#02ccfe", "#000000")

mycolors <- c(mycolors, 
			unclass(wes_palette(15, name="Zissou1", type="continuous")),
			unclass(wes_palette(15, name="Rushmore1", type="continuous")),
			unclass(wes_palette(15, name="FantasticFox1", type="continuous")))

mutation_colors <- c("#8effc1" = "Low",
					 "orange" = "Medium", 
					 "red" = "High")


# Functions for ordering x-axis (from stackoverflow)
reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
  new_x <- paste(x, within, sep = sep)
  stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
  reg <- paste0(sep, ".+$")
  ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}

# Make unique primer name column
primer_result <- primer_result %>%
  mutate(primer = paste0(ori_primer, "_", type))
primer_mutations <- primer_mutations %>%
  mutate(primer = paste0(ori_primer, "_", type))


## Grouped primer mutations
# plot_mutations_primers_take2 <- function(mydf) {
#   plotdf <- mydf %>%
# 	select(primer, primer_set, type, mutations, alert, count_genomes_with_mutation, count_genomes_amplified) %>%
# 	mutate(percentage = round(count_genomes_with_mutation/count_genomes_amplified*100, 2),
# 		   mylabel=paste0(count_genomes_with_mutation, "\n","(", percentage, "%)")) %>%
# 	rename(count = count_genomes_with_mutation) 
  
#   myplot <- plotdf %>%
# 	ggplot(aes(x=reorder_within(mutations, count, primer), y=count)) +
# 	geom_bar(aes(fill=alert, color=primer_set), stat="identity", position="stack") +
# 	geom_text(aes(label=mylabel, group=primer_set), position = position_stack(vjust=0.5), size=2.5) +
# 	facet_wrap(~primer, scales="free") +
# 	scale_fill_manual(values=names(mutation_colors), name="Alert status") +
# 	theme_bw() +
# 	scale_x_reordered() + 
# 	scale_y_continuous(expand=c(0.08, 0.05)) +
# 	ylab("Number of genomes") +
# 	xlab("Mismatches in primer binding sites") +
# 	theme(axis.text = element_text(size=7),
# 		  axis.title = element_text(size=8.5)) +
# 	coord_flip()
#   myplot
# }

plot_mutations_primers <- function(mydf) {
  plotdf <- mydf %>%
	select(-genomes_with_mutation, -primer_sequence, -binding_site) %>%
	group_by(primer, mutations, alert) %>%
	summarize(count_genomes_with_mutation=sum(count_genomes_with_mutation),
			  count_genomes_amplified = mean(count_genomes_amplified),
			  perc_genomes_amplified = mean(perc_genomes_amplified),
			  perc_with_mutation_amplified = sum(perc_with_mutation_amplified),
			  perc_with_mutation_total = sum(perc_with_mutation_total)) %>%
	mutate(percentage = round(count_genomes_with_mutation/count_genomes_amplified*100, 2),
		   mylabel=paste0(count_genomes_with_mutation, "\n","(", percentage, "%)")) %>%
	rename(count = count_genomes_with_mutation) 
  
	myplot <- plotdf %>%
	#group_by(primer, mutations, alert, count_genomes_amplified, )
		ggplot(aes(x=reorder_within(mutations, count, primer), y=count)) +
		geom_bar(aes(fill=alert), stat="identity") +
		geom_text(aes(label=mylabel), position = position_stack(vjust=0.5), size=2.5) +
		facet_wrap(~primer, scales="free") +
		scale_fill_manual(values=names(mutation_colors), name="Alert status") +
		theme_bw() +
		scale_x_reordered() + 
		scale_y_continuous(expand=c(0.08, 0.05)) +
		ylab("Number of genomes") +
		xlab("Mismatches in primer binding sites") +
		theme(axis.text = element_text(size=7),
			axis.title = element_text(size=8.5)) +
		coord_flip()
  myplot
}

primerplot_grouped <- plot_mutations_primers(primer_result)


## Individual primer mutations
# plot_individual_primer_mutations_take2 <- function(mydf) {
#   plotdf <- mydf %>%
# 	select(-genomes_with_mutation) %>%
# 	left_join(distinct(select(primer_result, primer, count_genomes_amplified)), by="primer") %>%
# 	group_by(position, primer) %>%
# 	mutate(count_wildtype = count_genomes_amplified - sum(count_genomes_with_mutation)) %>%
# 	select(-ori_primer, -type, -count_genomes_amplified) %>% 
# 	pivot_longer(c(count_genomes_with_mutation, count_wildtype), names_to="category", values_to="count") %>% 
# 	mutate(mutation=case_when(category=="count_wildtype" ~ "WT", TRUE ~ mutation)) %>%
# 	distinct() %>%
# 	ungroup() %>%
# 	group_by(primer, position) %>%
# 	mutate(percentage = count/sum(count)*100)
#   plotdf$position <- factor(plotdf$position, levels=unique(plotdf$position))

#   myplot <- plotdf %>%
# 	ggplot(aes(x=position, y=count)) +
# 	geom_bar(aes(fill=mutation), position="stack", stat="identity") +
# 	geom_text(aes(label=ifelse(percentage >=3, mutation, NA_character_), group=mutation),
# 					position = position_stack(vjust=0.5), max.overlaps=5, size=2.5, point.size=NA) +
# 	facet_wrap(.~primer, scales="free_x") + 
# 	scale_y_continuous(expand=c(0.08, 0.05)) +
# 	scale_fill_manual(values=mycolors) +
# 	theme_bw() +
# 	ylab("Number of genomes") +
# 	xlab("Position in primer sequence") +
# 	theme(legend.position = "blank", 
# 		  axis.text = element_text(size=8),
# 		  axis.title = element_text(size=8.5)) 
#   myplot
# }


plot_individual_primer_mutations <- function(mydf) {
	df_mut <- mydf %>%
		select(-genomes_with_mutation, -type, -primer_set) %>%
		group_by(ori_primer, primer, position, mutation) %>%
		summarise(count_genomes_with_mutation = sum(count_genomes_with_mutation)) 
	
	# Make df of wild types
	df_wt <- df_mut %>% 
		select(ori_primer, primer, position, count_genomes_with_mutation) %>%
		group_by(ori_primer, primer, position) %>%
		summarize(count_genomes_with_mutation = sum(count_genomes_with_mutation)) %>%
		left_join(distinct(select(primer_result, primer, count_genomes_amplified)), by="primer") %>%
		mutate(mutation="WT", count_genomes_wt = count_genomes_amplified - count_genomes_with_mutation) %>%
		relocate(mutation, .after=position) %>%
		select(-count_genomes_with_mutation, -count_genomes_amplified) %>%
		rename(count_genomes_with_mutation = count_genomes_wt) %>%
		distinct()
	  
	merged <- bind_rows(df_mut, df_wt) %>% 
		arrange(position) %>%
		group_by(primer, position) %>%
		mutate(percentage = count_genomes_with_mutation/sum(count_genomes_with_mutation)*100) 
	merged$position <- factor(merged$position, levels=unique(merged$position))
	
	myplot <- merged %>%
		ggplot(aes(x=position, y=count_genomes_with_mutation)) +
		geom_bar(aes(fill=mutation), position="stack", stat="identity") +
		geom_text(aes(label=ifelse(percentage >=3, mutation, NA_character_), group=mutation),
				max.overlaps=5, size=2.5, position = position_stack(vjust=0.5), point.size=NA) +
		facet_wrap(.~primer, scales="free_x") + 
		scale_y_continuous(expand=c(0.08, 0.05)) +
		scale_fill_manual(values=mycolors) +
		theme_bw() +
		ylab("Number of genomes") +
		xlab("Position in primer sequence") +
		theme(legend.position = "blank", 
			axis.text = element_text(size=8),
			axis.title = element_text(size=8.5)) 
	myplot
}

primerplot_ind <- plot_individual_primer_mutations(primer_mutations)


## Grouped probe mutations
# plot_mutations_probes_take2 <- function(mydf) {
#   plotdf <- mydf %>%
# 	select(ori_primer, ori_probe, mutations, alert, count_genomes_with_mutation, count_genomes_present) %>%
# 	#select(ori_probe, mutations, alert, count_genomes_with_mutation, count_genomes_present) %>%
# 	mutate(percentage = round(count_genomes_with_mutation/count_genomes_present*100, 2),
# 		   mylabel=paste0(count_genomes_with_mutation, "\n","(", percentage, "%)")) %>%
# 	rename(count = count_genomes_with_mutation) 
  
#   myplot <- plotdf %>%
# 	ggplot(aes(x=reorder_within(mutations, count, ori_primer), y=count)) +
# 	#ggplot(aes(x=reorder_within(mutations, count, ori_probe), y=count)) +
# 	geom_bar(aes(fill=alert, color=ori_probe), stat="identity", position="stack") +
# 	geom_text(aes(label=mylabel, group=ori_probe), position = position_stack(vjust=0.5), size=2.5) +
# 	facet_wrap(~ori_primer, scales="free") +
# 	scale_fill_manual(values=names(mutation_colors), name="Alert status") +
# 	theme_bw() +
# 	scale_x_reordered() + 
# 	scale_y_continuous(expand=c(0.08, 0.05)) +
# 	ylab("Number of genomes") +
# 	xlab("Mismatches in probe binding sites") +
# 	theme(axis.text = element_text(size=7),
# 		  axis.title = element_text(size=8.5)) +
# 	coord_flip()
#   myplot
# }

plot_mutations_probes <- function(mydf) {
  plotdf <- mydf %>%
    select(-genomes_with_mutation, -probe_sequence, -binding_site) %>%
    group_by(ori_primer, mutations, alert) %>%
    summarize(count_genomes_with_mutation = sum(count_genomes_with_mutation),
              count_genomes_present = mean(count_genomes_present),
              perc_genomes_present = mean(perc_genomes_present),
              perc_with_mutation_present = sum(perc_with_mutation_present),
              perc_with_mutation_total = sum(perc_with_mutation_total)) %>%
    mutate(percentage = round(count_genomes_with_mutation/count_genomes_present*100, 2),
           mylabel=paste0(count_genomes_with_mutation, "\n","(", percentage, "%)")) %>%
    rename(count = count_genomes_with_mutation) 
  
  myplot <- plotdf %>%
    ggplot(aes(x=reorder_within(mutations, count, ori_primer), y=count)) +
    #ggplot(aes(x=reorder_within(mutations, count, ori_probe), y=count)) +
    geom_bar(aes(fill=alert), stat="identity", position="stack") +
    geom_text(aes(label=mylabel), position = position_stack(vjust=0.5), size=2.5) +
    facet_wrap(~ori_primer, scales="free") +
    scale_fill_manual(values=names(mutation_colors), name="Alert status") +
    theme_bw() +
    scale_x_reordered() + 
    scale_y_continuous(expand=c(0.08, 0.05)) +
    ylab("Number of genomes") +
    xlab("Mismatches in probe binding sites") +
    theme(axis.text = element_text(size=7),
          axis.title = element_text(size=8.5)) +
    coord_flip()
  myplot
}

probeplot_grouped <- plot_mutations_probes(probe_result)
  
  

## Individual probe mutations
# plot_individual_probe_mutations_take2 <- function(mydf) {
#   plotdf <- mydf %>%
# 	select(-genomes_with_mutation) %>%
# 	left_join(distinct(select(probe_result, ori_probe, ori_primer, count_genomes_present)), by="ori_probe") %>%
# 	group_by(position, ori_primer) %>%
# 	#group_by(position, ori_probe) %>%
# 	mutate(count_wildtype = count_genomes_present - sum(count_genomes_with_mutation)) %>%
# 	select(-count_genomes_present) %>% 
# 	pivot_longer(c(count_genomes_with_mutation, count_wildtype), names_to="category", values_to="count") %>% 
# 	mutate(mutation=case_when(category=="count_wildtype" ~ "WT", TRUE ~ mutation)) %>%
# 	distinct() %>%
# 	ungroup() %>%
# 	group_by(ori_primer, position) %>%
# 	mutate(percentage = count/sum(count)*100)
#   plotdf$position <- factor(plotdf$position, levels=sort(unique(plotdf$position)))
  
#   myplot <- plotdf %>%
# 	ggplot(aes(x=position, y=count)) +
# 	geom_bar(aes(fill=mutation), position="stack", stat="identity") +
# 	geom_text(aes(label=ifelse(percentage >=3, mutation, NA_character_), group=mutation),
# 					position = position_stack(vjust=0.5), max.overlaps=5, size=2.5, point.size=NA) +
# 	facet_wrap(.~ori_primer, scales="free_x") + 
# 	scale_y_continuous(expand=c(0.01, 0.01)) +
# 	scale_fill_manual(values=mycolors) +
# 	theme_bw() +
# 	ylab("Number of genomes") +
# 	xlab("Position in probe sequence") +
# 	theme(legend.position = "blank", 
# 		  axis.text = element_text(size=8),
# 		  axis.title = element_text(size=8.5)) 
#   myplot
# }


plot_individual_probe_mutations <- function(mydf) {
	byprobe <- probe_grouped %>%
		group_by(ori_primer) %>%
		summarize(count_genomes_present = n_distinct(genome))
	
	df_mut <- mydf %>%
		select(-genomes_with_mutation, -type) %>%
		group_by(ori_primer, position, mutation) %>%
		summarise(count_genomes_with_mutation = sum(count_genomes_with_mutation))

	df_wt <- df_mut %>%
		group_by(ori_primer, position) %>%
		summarize(count_genomes_with_mutation = sum(count_genomes_with_mutation)) %>%
		left_join(distinct(select(byprobe, ori_primer, count_genomes_present)), by="ori_primer") %>%
		mutate(mutation="WT", count_genomes_wt = count_genomes_present - count_genomes_with_mutation) %>%
		relocate(mutation, .after=position) %>%
		select(-count_genomes_with_mutation, -count_genomes_present) %>%
		dplyr::rename(count_genomes_with_mutation = count_genomes_wt) %>%
		distinct()

	merged <- bind_rows(df_mut, df_wt) %>% 
		arrange(position) %>%
		group_by(ori_primer, position) %>%
		mutate(percentage = count_genomes_with_mutation/sum(count_genomes_with_mutation)*100) 
	merged$position <- factor(merged$position, levels=unique(merged$position))

	myplot <- merged %>%
		ggplot(aes(x=position, y=count_genomes_with_mutation)) +
		geom_bar(aes(fill=mutation), position="stack", stat="identity") +
		geom_text(aes(label=ifelse(percentage >=3, mutation, NA_character_), group=mutation),
					max.overlaps=5, size=2.5, position = position_stack(vjust=0.5), 
					point.size=NA) +
		facet_wrap(.~ori_primer, scales="free_x") + 
		scale_y_continuous(expand=c(0.01, 0.01)) +
		scale_fill_manual(values=mycolors) +
		theme_bw() +
		ylab("Number of genomes") +
		xlab("Position in probe sequence") +
		theme(legend.position = "blank", 
				axis.text = element_text(size=8),
				axis.title = element_text(size=8.5)) 
  myplot
}


probeplot_ind <- plot_individual_probe_mutations(probe_mutations)


plotprimer <- plot_grid(primerplot_grouped, primerplot_ind, nrow=2, 
			rel_heights=c(0.65, 0.35), labels=c("a", "b"))

plotprobe <- plot_grid(probeplot_grouped, probeplot_ind, nrow=2, 
			rel_heights=c(0.6, 0.4), labels=c("a", "b"))



# Alluvial plot of genome status
primer_status <- primer_status %>%
  mutate(key = paste0(genome, "_", ori_primer)) %>%
  rename(status_primer = status) 
probe_status <- probe_status %>%
  mutate(key = paste0(genome, "_", ori_primer)) %>%
  select(-ori_primer, -genome) %>%
  rename(status_probe = status)
status <- primer_status %>%
  left_join(probe_status, by="key") %>%
  select(-ori_probe) %>%
  group_by(status_primer, status_probe, ori_primer) %>%
  summarise(count=n())


alluvialplot <- ggplot(status, aes(axis1 = status_primer, axis2 = status_probe, y = count)) +
  geom_alluvium(aes(fill = status_primer)) + 
  geom_stratum(alpha=0.25) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3, color="black") +
  scale_x_discrete(limits = c("Primer status", "Probe status"), expand = c(0.15, 0.15)) +
  facet_wrap(~ori_primer) +
  theme_void() +
  theme(
	axis.text.x = element_text(angle = 0), 
	axis.title.x = element_text(size = 8),  
	axis.title.y = element_blank(),  
	strip.text = element_text(size = 8),  
	panel.grid = element_blank(),  
	legend.position = "none") 

ggsave(primeroutput, plotprimer, width=12, height=12)
ggsave(probeoutput, plotprobe, width=12, height=10)

ggsave(primer_combo_plot, primerplot_grouped, width=12, height=9)
ggsave(primer_mutation_plot, primerplot_ind, width=12, height=5)
ggsave(probe_combo_plot, probeplot_grouped, width=12, height=9)
ggsave(probe_mutation_plot, probeplot_ind, width=12, height=5)
ggsave(alluvial_plot, alluvialplot, width=4, height=7, bg="white")