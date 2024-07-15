#!/usr/bin/env Rscript

library(ggtree)
library(data.table)
library(ape)
library(dplyr)
library(tidyr)
library(stringr)
library(phytools)
library(RColorBrewer)
library(Polychrome)
library(ggplot2)
library(ggnewscale)
library(countrycode)


#$1: primer info file
#$2: genome list file
#$3: bed file
#$4: joined stats output file
#$5: target metadata file
#$6: full statistics output file
#$7: per primer statistics output file
#$8: per genome statistics output file
#$9: barplot output file
#$10: treeplot output file

args <- commandArgs(trailingOnly=TRUE)

primers <- fread(args[1], header=F, sep="\t")
target_genomes <- readLines(args[2])
bed <- fread(args[3], header=F, sep="\t")
tree <- read.tree(args[4])
target_metadata <- fread(args[5], header=T, sep="\t")
accession <- fread(args[6], header=F, sep="\t")
primers_ori <- fread(args[7], header=F, sep="\t")

## Output file names
full_stats <- args[8]
primer_stats <- args[9]
genome_stats <- args[10]
barplot <- args[11]
treeplot <- args[12]


genomecount <- length(target_genomes)

if (!is.null(accession)) {
	names(accession) <- c("genome", "assembly_id")
}

names(primers) <- c("genome", "start", "end", "primerset", "product_size", "fwd_primer", "rev_primer")
names(primers_ori) <- c("primerset_ori", "primerset", "fwd_primer", "rev_primer", "probe")
names(bed) <- c("genome", "start", "end", "primerset", "score", "strand")

## Lookup assembly id for genome accessions
primers <- primers %>% 
	left_join(accession, by="genome") %>%
	select(-genome) %>%
	rename(genome=assembly_id) %>%
	relocate(genome)

bed <- bed %>%
	left_join(accession, by="genome") %>%
	select(-genome) %>%
	rename(genome=assembly_id) %>%
	relocate(genome)

## Get lengths of primers
primers <- primers %>%
	mutate(sum_primer_lengths = nchar(fwd_primer) + nchar(rev_primer))
primers$product_size <- gsub("bp", "", primers$product_size)
primers <- primers %>%
	mutate(joincol = paste(genome, primerset, start, end, sep="_"))

## Left join primer info to bed file
## Need score and strand info from bed file
bed$start <- bed$start+1

bed <- bed %>%
	mutate(joincol = paste(genome, primerset, start, end, sep="_"))
primers <- primers %>%
	left_join(select(bed, joincol, score, strand), by="joincol") %>%
	mutate(sum_mismatches = round(sum_primer_lengths-(sum_primer_lengths*score/1000)),
		   mismatch_percentage = sum_mismatches/sum_primer_lengths*100) %>%
	left_join(select(primers_ori, primerset, primerset_ori), by="primerset")

fwrite(bed, file=full_stats, sep="\t", quote=F, col.names=T, row.names=F)

## Process target metadata
## if using nt
if ("Geographic Location" %in% names(target_metadata)) {
	target_metadata <- target_metadata %>%
		select(Accession, `Geographic Location`, `Geographic Region`, `Isolate Collection date`, `Host Name`) %>%
		rename(genome = Accession, 
				country = `Geographic Location`, 
				continent = `Geographic Region`, 
				year = `Isolate Collection date`,
				host = `Host Name`) %>%
		mutate(across(everything(), ~replace_na(.x, "NA"))) %>%
		# replace Viet Nam with Vietnam
		mutate(country=case_when(country=="Viet Nam" ~ "Vietnam", TRUE ~ country)) %>%
		separate(country, into=c("country", "discard"), sep=":", remove=T, extra="merge") %>%
		select(-discard) %>%
		# extract first 4 characters in year column
		mutate(year=str_sub(year, start=1, end=4)) %>%
		filter(genome %in% tree$tip.label) %>%
		mutate(across(everything(), ~ifelse(.=="", NA, as.character(.)))) %>%
		rename(id=genome) %>%
		mutate(genome=id) %>%
		distinct()
## if using assembly
} else {
	target_metadata <- target_metadata %>%
		filter(`Assembly Accession` %in% tree$tip.label) %>%
		rename(genome = `Assembly Accession`,
				country = `Assembly BioSample Location name`,
				host = `Assembly BioSample Host`,
				year = `Assembly BioSample Collection date`,
				label = `Organism Infraspecific Names Strain`) %>%
		separate(country, into=c("country", "discard"), sep=":", remove=T, extra="merge") %>%
		mutate(country=case_when(country=="Viet Nam" ~ "Vietnam", TRUE ~ country)) %>%
		select(-discard) %>%
		mutate(year = str_sub(year, start=1, end=4)) %>%
		
		mutate(across(everything(), ~ifelse(.=="", NA, as.character(.)))) %>%
		rename(id=genome) %>%
		mutate(genome=label) %>% 
		select(-label)
	target_metadata$continent <- countrycode(target_metadata$country, origin="country.name", destination="continent")
	target_metadata <- target_metadata %>%
		relocate(continent, .after=country) %>%
		distinct()
}


# Generate per primerset summary statistics
summary <- primers %>%
	group_by(primerset_ori) %>%
	summarize(`Perfect match` = n_distinct(genome[sum_mismatches==0])/genomecount, 
			  `≤1 mismatch` = n_distinct(genome[sum_mismatches<=1])/genomecount,
			  `≤2 mismatches` = n_distinct(genome[sum_mismatches<=2])/genomecount,
			  `≤3 mismatches` = n_distinct(genome[sum_mismatches<=3])/genomecount,
			  `≤4 mismatches` = n_distinct(genome[sum_mismatches<=4])/genomecount,
			  `≤5 mismatches` = n_distinct(genome[sum_mismatches<=5])/genomecount,
			  `≤10 mismatches` = n_distinct(genome[sum_mismatches<=10])/genomecount)


fwrite(summary, file=primer_stats, sep="\t", quote=F, col.names=T, row.names=F)

summary_reshape <- summary %>%
	tidyr::pivot_longer(-primerset_ori, names_to="category", values_to="prevalence")

## Generate per genome statistics
summary_genome_mismatches <- primers %>%
	# exclude hits with more than 10 mismatches
	filter(sum_mismatches <= 10) %>%
	select(genome, primerset_ori, sum_mismatches) %>%
	group_by(genome, primerset_ori) %>%
	summarize(sum_mismatches=min(sum_mismatches)) %>%
	mutate(sum_mismatches=case_when(sum_mismatches>5 ~ 6, TRUE ~ sum_mismatches)) %>%
	ungroup() %>%
	complete(genome, primerset_ori, fill=list(sum_mismatches=NA)) %>%
	tidyr::pivot_wider(names_from="primerset_ori", 
	values_from="sum_mismatches")

## Get list of genomes not present in summary file
absent <- setdiff(target_genomes, summary_genome_mismatches$genome)

## Create tmp data frame with absent genomes
tmp_df <- setNames(data.frame(matrix(ncol = ncol(summary_genome_mismatches), nrow = length(absent))), names(summary_genome_mismatches))
tmp_df$genome <- absent

## Merge tmp df into summary df
summary_genome_mismatches <- bind_rows(summary_genome_mismatches, tmp_df)

## Df for generating heatmap
summary_genome_forheatmap <- primers %>%
	select(genome, primerset_ori, sum_mismatches) %>%
	complete(genome, primerset_ori, fill=list(sum_mismatches=NA)) 

fwrite(summary_genome_mismatches, file=genome_stats, sep="\t", quote=F, col.names=T, row.names=F)


## Plot figures by primerset, faceted by mismatch count
summary_reshape$category <- factor(summary_reshape$category, 
	levels=c("Perfect match", "≤1 mismatch", "≤2 mismatches", "≤3 mismatches", "≤4 mismatches", "≤5 mismatches", "≤10 mismatches"))

primers.count <- length(unique(summary_reshape$primerset_ori))
primers.colors <- unname(sample(palette36.colors()[names(palette36.colors())!="Purplish_White"], size=primers.count))

myplot <- summary_reshape %>%
	ggplot(aes(x=category, y=prevalence)) +
	geom_bar(stat="identity", position="dodge", aes(fill=category)) +
	facet_wrap(.~primerset_ori, nrow=1) +
	ylab("Percentage of genomes") +
	scale_fill_manual(values=primers.colors, name="Primer set") +
	scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
	theme(axis.text.x=element_blank(), 
		  axis.title.x=element_blank(),
		  axis.ticks.x=element_blank())

ggsave(barplot, myplot, width=8, height=6, dpi = 300, device="pdf")

## Plot per genome heatmap
# myheatmap <- summary_genome_forheatmap %>%
# 	ggplot(aes(x=primerset, y=genome)) +
# 	geom_tile(aes(fill=sum_mismatches), width=0.5, height=0.5, size=0.8, color="#909090") +
# 	scale_fill_continuous_sequential(palette="Emrld", na.value="grey85") +
# 	#scale_fill_gradient(low="#188300", high="#d1fbd8", na.value = "grey90") +
# 	labs(x="Primer set", y="Genome", fill="Sum primer mismatches") +
# 	scale_y_discrete(expand=c(0,0)) +
# 	theme_minimal() +
# 	theme(axis.text.x=element_text(angle=45, size=8),
# 		 axis.text.y=element_text(size=5.5))

# ggsave(heatmap, myheatmap, width=5, height=20, dpi = 300, device="pdf")


##############################################################################################################################
## Plot tree
## First root tree at midpoint if unrooted

if (!is.rooted(tree)) {
	tree.mp <- phytools::midpoint(tree)
} else {
	tree.mp <- tree
}

## Mismatches data for heatmap
heatmap_data <- summary_genome_mismatches
rownames(heatmap_data) <- heatmap_data$genome
rn <- rownames(heatmap_data)
heatmap_data <- heatmap_data[-1]

heatmap_data <- as.data.frame(sapply(heatmap_data, as.character))
rownames(heatmap_data) <- rn
heatmap_data <- heatmap_data %>% 
	filter(rownames(.) %in% tree$tip.label)

## Other sample metadata for annotation
target_metadata2 <- target_metadata
rownames(target_metadata2) <- target_metadata2$id
rn2 <- rownames(target_metadata2)
target_metadata2 <- target_metadata2[-1]
target_metadata2 <- as.data.frame(sapply(target_metadata2, as.character))
rownames(target_metadata2) <- rn2

## Define colors for sample metadata
country.count <- length(unique(target_metadata$country))
continent.count <- length(unique(target_metadata$continent))
year.count <- length(unique(target_metadata$year))
host.count <- length(unique(target_metadata$host))
#cols <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

country.colours <- unname(sample(palette36.colors()[names(palette36.colors())!="Purplish_White"], size=country.count))
continent.colours <- unname(sample(kelly.colors()[names(kelly.colors())!="white"], size=continent.count))
host.colours <- unname(sample(palette36.colors()[names(palette36.colors())!="Purplish_White"], size=country.count))
year.colours <- colorRampPalette(brewer.pal(8, "Reds"))(year.count)
mismatch.colours <- c("#74c365", "#006b3e", "#ffe733", "#FF8C01", "#ff0800", "#7c0a02", "#300401")

## Function to make legends smaller to fit in page
addSmallLegend <- function(myplot, pointSize=11, textSize=7.5, spaceLegend=0.8) {
  myplot +
	guides(shape = guide_legend(override.aes = list(size = pointSize)),
		   color = guide_legend(override.aes = list(size = pointSize))) +
	theme(legend.title = element_text(size = 13), 
		  legend.text  = element_text(size = textSize),
		  legend.key.size = unit(spaceLegend, "lines"))
}

# Plot tree with continent as colored tips
p <- ggtree(tree.mp, layout = "circular", size=0.2) %<+% target_metadata +
  geom_tippoint(aes(color=continent), size=0.5) + 
  scale_color_manual(values=continent.colours, name="Continent") +
  geom_tiplab(aes(label=genome), 
			  align=F, linetype=NA, size=1, 
			  offset=0.13, hjust=0.5) 

# Plot country as heatmap
p1 <- gheatmap(p, (target_metadata2 %>% select(country)), 
	offset = 0.2, 
	colnames_position="top", 
	colnames_angle=90, 
	colnames_offset_y = 0.2,
	width=0.15, hjust=0, 
	font.size=1.8, 
	color="#909090") +
  scale_fill_manual(values=country.colours, 
					na.value="#d3d3d3", name="Country")
p1tmp <- addSmallLegend(p1) + new_scale_fill()

# Plot year as heatmap
p2 <- gheatmap(p1tmp, (target_metadata2 %>% select(year)), 
	offset = 0.33, 
	colnames_position="top", 
	colnames_angle=90, 
	colnames_offset_y = 0.33, 
	width=0.15, hjust=0, 
	font.size=1.8, 
	color="#909090") +
  scale_fill_manual(values=year.colours, 
					na.value="#d3d3d3", name="Year")
p2tmp <- addSmallLegend(p2) + new_scale_fill()

# Plot host as heatmap
p3 <- gheatmap(p2tmp, (target_metadata2 %>% select(host)), 
	offset = 0.46, 
	colnames_position="top", 
	colnames_angle=90, 
	colnames_offset_y = 0.46, 
	width=0.15, hjust=0, 
	font.size=1.8, 
	color="#909090") +
  scale_fill_manual(values=host.colours, 
					na.value="#d3d3d3", name="Host")
p3tmp <- addSmallLegend(p3) + new_scale_fill()

# Plot mismatches as heatmap
width_mismatches <- 0.1 * nrow(summary)
p4 <- gheatmap(p3tmp, heatmap_data, 
	offset = 0.59, 
	colnames_position="top", 
	colnames_angle=90, 
	colnames_offset_y = 0.59, 
	width=width_mismatches, hjust=0, 
	font.size=1.8, 
	color="#909090") +
  scale_fill_manual(values=mismatch.colours, 
  					na.value="#d3d3d3", name="Overall primer mismatch count")
p4 <- addSmallLegend(p4)

ggsave(treeplot, plot=p4, device="pdf", height=20, width=20)
