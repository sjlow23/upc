#!/usr/bin/env Rscript

library(data.table)
library(igraph)
library(ape)
library(phangorn)

args <- commandArgs(trailingOnly=TRUE)

data <- fread(args[1], sep="\t", header=F)
names(data) <- c("g1", "g2", "dist")

outtree <- args[2]

g <- graph_from_data_frame(data, directed=T)
g2 <- as.undirected(g, mode="collapse", edge.attr.comb="mean")

adj_matrix <- as.matrix(as_adjacency_matrix(g2, attr="dist", type="both", sparse = FALSE))
dist_matrix <- as.dist(adj_matrix)
tree <- nj(dist_matrix)
tree.mp <- midpoint(tree)
write.tree(tree.mp, file=outtree)


