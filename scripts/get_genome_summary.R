#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)

args <- commandArgs(trailingOnly=TRUE)


mismatch <- fread(args[1], sep="\t", header=T)
missing <- fread(args[2], sep="\t", header=T)
output <- args[4]
genomelist <- gsub(".fna", "", readLines(args[5]))
type <- args[6]

if (type == "probes") {
    query <- unique(fread(args[3], sep="\t", header=F)$V2)
    probes_lookup <- fread(args[3], sep="\t", header=F)
    names(probes_lookup) <- c("primer", "probe", "probe_seq")
    # Get genomes with mismatches
    mismatchdf <- mismatch %>%
        select(genome, ori_probe) %>%
        distinct() %>%
        rename(probe = ori_probe) %>%
        left_join(select(probes_lookup, probe, primer), by="probe") %>%
        mutate(status = "with mismatch")

    # Get missing genomes
    missingdf <- missing %>%
        filter(count_genomes_missing != 0) %>%
        select(ori_probe, genomes_missing) %>%
        separate_rows(genomes_missing, sep=", ") %>%
        mutate(status = "not present") %>%
        rename(genome = genomes_missing, probe = ori_probe) %>%
        left_join(select(probes_lookup, probe, primer), by="probe") %>%
        relocate(genome)

    # For each probe set, identify genomes with perfect matches
    get_perfectmatches <- function(myprobe) {
        mismatch <- mismatchdf %>%
            filter(probe == myprobe) %>%
            pull(genome)
        missing <- missingdf %>%
            filter(probe == myprobe) %>%
            pull(genome)
        perfect <- setdiff(genomelist, c(mismatch, missing))
        perfectdf <- data.frame(genome=perfect, probe=myprobe, status="perfect match")
        perfectdf <- perfectdf %>%
            left_join(select(probes_lookup, probe, primer), by="probe") %>%
            relocate(status, .before="primer")
    }
}


if (type == "primers") {
    query <- unique(fread(args[3], sep="\t", header=F)$V1)
    # Get genomes with mismatches
    mismatchdf <- mismatch %>%
        select(genome, ori_primer) %>%
        distinct() %>%
        rename(primer = ori_primer) %>%
        mutate(status = "with mismatch")

    # Get missing genomes
    missingdf <- missing %>%
        filter(count_genomes_missing != 0) %>%
        select(ori_primer, genomes_missing) %>%
        separate_rows(genomes_missing, sep=", ") %>%
        mutate(status = "not amplified") %>%
        rename(genome = genomes_missing, primer = ori_primer) %>%
        relocate(genome)

    # For each primer set, identify genomes with perfect matches
    get_perfectmatches <- function(myprimer) {
        mismatch <- mismatchdf %>%
            filter(primer == myprimer) %>%
            pull(genome)
        missing <- missingdf %>%
            filter(primer == myprimer) %>%
            pull(genome)
        perfect <- setdiff(genomelist, c(mismatch, missing))
        perfectdf <- data.frame(genome=perfect, primer=myprimer, status="perfect match")
    }
}


combined_perfect <- rbindlist(lapply(query, get_perfectmatches))

merged_res <- bind_rows(mismatchdf, missingdf, combined_perfect)
merged_res <- merged_res %>%
    arrange(.[[1]], .[[2]]) %>%
    distinct()

fwrite(merged_res, output, col.names=T, row.names=F, sep="\t", quote=F)