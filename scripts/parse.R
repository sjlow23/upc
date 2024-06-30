#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(colorspace)
library(RColorBrewer)
library(scales)
library(ape)
library(ggtree)
library(tidytree)
library(tibble)
library(ggnewscale)
library(stringr)


#$1: primer info file
#$2: genome list file
#$3: bed file
#$4: joined stats output file
#$5: per primer statistics output file
#$6: per genome statistics output file
#$7: barplot output file
#$8: heatmap output file

args <- commandArgs(trailingOnly=TRUE)

primers <- fread(args[1], header=F, sep="\t")
target_genomes <- readLines(args[2])
bed <- fread(args[3], header=F, sep="\t")
full_stats <- args[4]
primer_stats <- args[5]
genome_stats <- args[6]
barplot <- args[7]
heatmap <- args[8]
tree <- read.tree(args[9])
treeplot <- args[10]
target_metadata <- fread(args[11], header=T, sep="\t")

genomecount <- length(target_genomes)

names(primers) <- c("genome", "start", "end", "primerset", "product_size", "fwd_primer", "rev_primer")
names(bed) <- c("genome", "start", "end", "primerset", "score", "strand")

## Get lengths of primers
primers <- primers %>%
	mutate(sum_primer_lengths = nchar(fwd_primer) + nchar(rev_primer))
primers$product_size <- gsub("bp", "", primers$product_size)

## Left join primer info to bed file
bed <- bed %>%
	select(-start, -end) %>%
	left_join(primers, by=c("genome", "primerset")) %>%
	mutate(sum_mismatches = round(sum_primer_lengths-(sum_primer_lengths*score/1000)),
		   mismatch_percentage = sum_mismatches/sum_primer_lengths*100)

fwrite(bed, file=full_stats, sep="\t", quote=F, col.names=T, row.names=F)

## Process target metadata
target_metadata <- target_metadata %>%
	select(Accession, `Geographic Location`, `Geographic Region`, `Isolate Collection date`) %>%
	rename(genome=Accession, 
			country=`Geographic Location`, 
			continent=`Geographic Region`, 
			year=`Isolate Collection date`) %>%
	mutate(across(everything(), ~replace_na(.x, "NA"))) %>%
	# replace Viet Nam with Vietnam
	mutate(country=case_when(country=="Viet Nam" ~ "Vietnam", TRUE ~ country)) %>%
	separate(country, into=c("country", "discard"), sep=":", remove=T) %>%
	select(-discard) %>%
	# extract first 4 characters in year column
	mutate(year=str_sub(year, start=1, end=4))
	

# Generate per primerset summary statistics
summary <- bed %>%
	group_by(primerset) %>%
	summarize(`Perfect match` = n_distinct(genome[sum_mismatches==0])/genomecount, 
			  `1 mismatch` = n_distinct(genome[sum_mismatches<=1])/genomecount,
			  `2 mismatches` = n_distinct(genome[sum_mismatches<=2])/genomecount,
			  `3 mismatches` = n_distinct(genome[sum_mismatches<=3])/genomecount,
			  `4 mismatches` = n_distinct(genome[sum_mismatches<=4])/genomecount,
			  `5 mismatches` = n_distinct(genome[sum_mismatches<=5])/genomecount)


fwrite(summary, file=primer_stats, sep="\t", quote=F, col.names=T, row.names=F)

summary_reshape <- summary %>%
	tidyr::pivot_longer(-primerset, names_to="category", values_to="prevalence")

## Generate per genome statistics
summary_genome_mismatches <- bed %>%
	select(genome, primerset, sum_mismatches) %>%
	complete(genome, primerset, fill=list(sum_mismatches=NA)) %>% 
	tidyr::pivot_wider(names_from="primerset", 
	values_from="sum_mismatches")

## Get list of genomes not present in summary file
absent <- setdiff(target_genomes, summary_genome_mismatches$genome)

## Create tmp data frame with absent genomes
tmp_df <- setNames(data.frame(matrix(ncol = ncol(summary_genome_mismatches), nrow = length(absent))), names(summary_genome_mismatches))
tmp_df$genome <- absent

## Merge tmp df into summary df
summary_genome_mismatches <- bind_rows(summary_genome_mismatches, tmp_df)

## Df for generating heatmap
summary_genome_forheatmap <- bed %>%
	select(genome, primerset, sum_mismatches) %>%
	complete(genome, primerset, fill=list(sum_mismatches=NA)) 

fwrite(summary_genome_mismatches, file=genome_stats, sep="\t", quote=F, col.names=T, row.names=F)


## Plot figures by primerset, faceted by mismatch count
summary_reshape$category <- factor(summary_reshape$category, 
	levels=c("Perfect match", "1 mismatch", "2 mismatches", "3 mismatches", "4 mismatches", "5 mismatches"))

myplot <- summary_reshape %>%
	ggplot(aes(x=category, y=prevalence)) +
	geom_bar(stat="identity", position="dodge", aes(fill=category)) +
	facet_wrap(.~primerset, nrow=1) +
	ylab("Percentage of genomes") +
	scale_fill_brewer(palette="Set2") +
	scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
	theme(axis.text.x=element_blank(), 
		  axis.title.x=element_blank(),
		  axis.ticks.x=element_blank())

ggsave(barplot, myplot, width=8, height=6, dpi = 300, device="pdf")

## Plot per genome heatmap
myheatmap <- summary_genome_forheatmap %>%
	ggplot(aes(x=primerset, y=genome)) +
	geom_tile(aes(fill=sum_mismatches), width=0.5, height=0.5, size=0.8) +
	scale_fill_continuous_sequential(palette="Emrld", na.value="grey85") +
	#scale_fill_gradient(low="#188300", high="#d1fbd8", na.value = "grey90") +
	labs(x="Primer set", y="Genome", fill="Sum primer mismatches") +
	scale_y_discrete(expand=c(0,0)) +
	theme_minimal() +
	theme(axis.text.x=element_text(angle=45, size=8),
		 axis.text.y=element_text(size=5.5))

ggsave(heatmap, myheatmap, width=5, height=20, dpi = 300, device="pdf")


## Plot tree
metadata <- summary_genome_mismatches

metadata <- metadata %>% 
	left_join(target_metadata, by="genome") %>%
	relocate(genome) %>%
	relocate(country, .after=genome)

# Define functions to plot tree with heatmap
# Credit to: http://dmnfarrell.github.io/r/ggtree-heatmaps

gettreedata <- function(tree, metadata) {
	d <- metadata[match(tree$tip.label, metadata$genome), ]
	d <- column_to_rownames(metadata, var="genome")
	d$label <- row.names(d)
	y <- full_join(as_tibble(tree), d, by='label')
	y <- as.treedata(y)
	return(y)
}

get_color_mapping <- function(data, col, cmap){
	labels <- (data[[col]])
	names <- levels(as.factor(labels))
	n <- length(names)
	if (n<10){
		colors <- suppressWarnings(c(brewer.pal(n, cmap)))[1:n]
	}
	else {
		colors <- colorRampPalette(brewer.pal(8, cmap))(n)
	}
	names(colors) = names
	return (colors)
}


ggplottree <- function(tree, meta, cols=NULL, colors=NULL, cmaps=NULL, layout="rectangular",
					   offset=10, tiplabel=FALSE, tipsize=3, tiplabelsize=5, tiplabelcol=NULL,
					   align=FALSE, tipoffset=0) {

	y <- gettreedata(tree, meta)
	if (layout == 'cladogram'){
		p <- ggtree(y, layout='c', branch.length='none')
	}
	else {
		p <- ggtree(y, layout=layout)
	}

	if (is.null(cols)) {
		if (tiplabel){
			p <- p + geom_tiplab(size=tiplabelsize,offset=tipoffset)
		}
		return (p)
	}
	col <- cols[1]
	if (!is.null(colors)) {
		#use predefined colors
		clrs <- colors
	}
	else {
		#calculate colors from cmap
		cmap <- cmaps[1]
		df <- meta[tree$tip.label,][col]
		clrs <- get_color_mapping(df, col, cmap)
	}
	#print (clrs)
	p <- p + new_scale_fill() +
			geom_tippoint(mapping=aes(fill=.data[[col]]),size=tipsize,shape=21,stroke=0) +
			scale_fill_manual(values=clrs, na.value="black")

	p2 <- p
	if (length(cols)>1){
		for (i in 2:length(cols)){
			col <- cols[i]
			if (length(cmaps)>=i){
				cmap <- cmaps[i]
			}
			else {
				cmap = 'Greys'
			}
			df <- meta[tree$tip.label,][col]
			type <- class(df[col,])
			print (type)
			p2 <- p2 + new_scale_fill()
			p2 <- gheatmap(p2, df, offset=i*offset, width=.08,
					  colnames_angle=0, colnames_offset_y = .05)
			if (type %in% c('numeric','integer')){
				p2 <- p2 + scale_fill_gradient(low='#F8F699',high='#06A958', na.value="white")
			}
			else {
				colors <- get_color_mapping(df, col, cmap)
				p2 <- p2 + scale_fill_manual(values=colors, name=col, na.value="white")
			}
		}
	}

	p2 <- p2 + theme_tree2(legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
						legend.position="left", plot.title = element_text(size=40))
			guides(color = guide_legend(override.aes = list(size=10)))
	if (tiplabel) {
		if (!is.null(tiplabelcol)) {
			p2 <- p2 + geom_tiplab(mapping=aes(label=.data[[tiplabelcol]]),
								size=tiplabelsize, align=align,offset=tipoffset)
		}
		else {
			p2 <- p2 + geom_tiplab(size=tiplabelsize, align=align,offset=tipoffset)
		}
	}
	return(p2)
}


treeplot <- ggplottree(tree, metadata, layout='rect', cols=colnames(metadata[-1]),
		cmaps=c('Set1','Set2','Blues'), tiplabel=TRUE, tipoffset=.1, tipsize=4, offset=.5)

ggsave(treeplot, finaltree, width=12, height=13))
