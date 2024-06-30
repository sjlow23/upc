library(ape)
library(ggtree)
library(ggplot2)
library(data.table)
library(tibble)

library(ggnewscale)
library(tidyverse)


#$1: tree file
#$2: metadata file

#$3: genome count (target or offtarget)
#$4: joined stats output file
#$5: per primer statistics output file
#$6: per genome statistics output file
#$7: barplot output file
#$8: heatmap output file

args <- commandArgs(trailingOnly=TRUE)

mytree <- read.tree(args[1])
metadata <- fread(args[2], header=T, sep="\t")
# metadata consist of genome in col1, primer names in col2 onwards (genomestats.tsv file)

# Keep only what's present in tree
#metadata <- metadata[match(mytree$tip.label, metadata$genome), ]

# Note: Not all genomes present in tree will be in metadata file
# Some genomes will not have any primers mapped, and will be lost from metadata file
# But they will still be in the tree
# Workaround: either prune tree to subset of common genomes with metadata file,
# or run genome step first, then infer tree
# can ggtree handle matrix with missing tree leaf data?


# Make sure tip label is first column in metadata file
metadata <- metadata %>% relocate(genome)

# Add rownames
#metadata2 <- column_to_rownames(metadata, var = "genome")

# Define functions to plot tree with heatmap
# Credit to: http://dmnfarrell.github.io/r/ggtree-heatmaps

gettreedata <- function(tree, metadata) {
    #get treedata object
    d <- metadata[match(mytree$tip.label, metadata$genome), ]
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

ggplottree <- function(tree, meta, cols=NULL, cmaps=NULL, layout="rectangular",
                       offset=10, tiplabel=FALSE, tipsize=3) {

    y <- gettreedata(tree, meta)
    p <- ggtree(y, layout=layout)
    if (is.null(cols)){
        return (p)
    }

    col <- cols[1]
    cmap <- cmaps[1]
    df<-meta[tree$tip.label,][col]
    colors <- get_color_mapping(df, col, cmap)

    #tip formatting
    p1 <- p + new_scale_fill() +
          geom_tippoint(mapping=aes(fill=.data[[col]]),size=tipsize,shape=21,stroke=0) +
          scale_fill_manual(values=colors, na.value="white")

    p2 <- p1
    if (length(cols)>1){
        for (i in 2:length(cols)){
            col <- cols[i]
            cmap <- cmaps[i]
            df <- meta[tree$tip.label,][col]
            type <- class(df[col,])
            p2 <- p2 + new_scale_fill()
            p2 <- gheatmap(p2, df, offset=i*offset, width=.08,
                      colnames_angle=0, colnames_offset_y = .05)
            #deal with continuous values
            if (type == 'numeric'){
                p2 <- p2 + scale_color_brewer(type="div", palette=cmap)
            }
            else {
                colors <- get_color_mapping(df, col, cmap)
                p2 <- p2 + scale_fill_manual(values=colors, name=col)
            }
        }
    }

    p2 <- p2 + theme_tree2(legend.text = element_text(size=20), legend.key.size = unit(1, 'cm'),
                        legend.position="left", plot.title = element_text(size=40))
            guides(color = guide_legend(override.aes = list(size=10)))

    return(p2)
}


treeplot <- ggplottree(mytree, metadata, cols=colnames(metadata[-1]),
           cmaps=c('Set1','Set2','Set3'), tipsize=8, offset=.5 ,layout='rect')


ggsave(snakemake@output[[3]], finaltree, width=12, height=13))









# Convert metadata dataframe to a list
# metadata2 <- metadata %>%
# 	pivot_longer(-genome, names_to="primerset", values_to="sum_mismatches")
# metadata_list <- as.list(metadata)

# # Create ggtree object from the tree
# mytree_g <- ggtree(mytree)


# ## Alternate view (aligning data based on tree structure)
# ## Visualize sum primer mismatches
# # Function to generate matrices for tree
# # select_meta_subset <- function(myfield) {
# # 	subsetmeta <- metadata %>%
# # 		dplyr::select(genome, myfield) %>% 
# # 		filter(genome %in% mytree$tip.label) %>%
# # 		tibble::remove_rownames() %>% 
# # 		tibble::column_to_rownames(var="genome")
# # 	subsetmeta
# # }

# # Make genomes row names for ggtree format
# metadata <- metadata %>% 
# 	column_to_rownames(var="genome")

# primersets <- unique(metadata2$primerset)

# #subset_all <- rbindlist(lapply(primersets, select_meta_subset))
# ggsave("annotated_tree.png", tree_annot, width = 8, height = 6, dpi = 300)

# library(ggtree)
# library(ggnewscale)
# library(tidyverse)
# library(data.table)


# #subsetmeta1 <- select_meta_subset("Country")
# #subsetmeta2 <- select_meta_subset("Year")
# #subsetmeta3 <- select_meta_subset("Host")

# # Function to plot alternative view of tree
# # Color tips- use other than mismatch color, maybe lineage or genotype, or number of primer sets with amplification
# plottree_alternative <- function() {
# 	treeplot2 <- ggtree(mytree) %<+% metadata2 + 
# 		geom_tippoint(aes(color=sum_mismatches), size=1.0)
# 	treeplot2
# }


# treeplot2 <- plottree_alternative()
# #gheatmap(treeplot2, subsetmeta, offset=1, width=0.5, colnames=F) 






# # Add matrices sequentially
# h1 <- gheatmap(treeplot2, p1, offset = 0, width = 0.05,
# 							 colnames_position = "top",
# 							 colnames_offset_y = 0.4,
# 							 font.size=2) +
# 	scale_fill_viridis(option = "D", name = "Sum mismatches") + coord_cartesian(clip = "off")

# h2 <- h1 + new_scale_fill()
# h3<- gheatmap(h2, subsetmeta2, offset = 0.01, width = 0.05,
# 							colnames_position = "top",
# 							colnames_offset_y = 0.5,
# 							font.size=2) +
# 	scale_fill_viridis_d(option = "D", name = "Year") + coord_cartesian(clip = "off")

# h4 <- h3 + new_scale_fill()  

# finaltree <- gheatmap(h4, subsetmeta3, offset = 0.02, width = 0.05,
# 				 colnames_position = "top", 
# 				 colnames_offset_y = 0.6,
# 				 font.size=2) +
# 	scale_fill_viridis_d(option = "D", name = "Host") + coord_cartesian(clip = "off")


# ggsave(snakemake@output[[3]], finaltree, width=12, height=13)
