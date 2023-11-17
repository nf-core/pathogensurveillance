#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- c("PIRATE.gene_families.ordered.tsv", "metadata_medium.csv", "10", "0.9", "0.5")
names(args) <- c("gene_families", "metadata", "min_core_genes", "min_core_samps", "min_core_refs")
args <- as.list(args)
gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = TRUE, sep = ',')
min_core_genes <- as.integer(args$min_core_genes)

# Infer number of samples and references
total_count <- ncol(gene_data) - 22
all_ids <- colnames(gene_data)[23:ncol(gene_data)]
sample_ids <- all_ids[all_ids %in% gsub(pattern = '[-.]+', replacement = '_', metadata$sample)] #TODO: change this once sample IDs are standardized
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Get minimum number of samples and references that is acceptable
min_core_samps <- ceiling(as.numeric(args$min_core_samps) * length(sample_ids))
min_core_refs <- ceiling(as.numeric(args$min_core_refs) * length(ref_ids))

# Remove samples until a core genome is found
find_worst_sample <- function(data) {
    
}

get_n_core_samples <- function(data) {
    rowSums(data[, sample_ids] != "")
}

get_n_core_refs <- function(data) {
    rowSums(data[, ref_ids] != "")
}

gene_data_subset <- gene_data
current_sample_ids <- sample_ids
current_ref_ids <- ref_ids
while (sum(get_n_core_samples(gene_data_subset) == length(current_sample_ids) & 
           get_n_core_refs(gene_data_subset) == length(current_ref_ids)) < min_core_genes) {
    
    # Find how many genes each sample is missing
    n_missing_samp <- colSums(gene_data_subset[, current_sample_ids] == "")
    n_missing_ref <- colSums(gene_data_subset[, current_ref_ids] == "")
    
    # Remove worst sample/ref if possible, preferring to remove references
    if (length(current_sample_ids) > min_core_samps && length(current_ref_ids) > min_core_refs) {
        worst_id <- names(which.min(c(n_missing_samp, n_missing_ref)))
    } else (length(current_ref_ids) > min_core_refs) {
        worst_id <- names(which.min(n_missing_ref))
    } else (length(current_sample_ids) > min_core_samps) {
        worst_id <- names(which.min(n_missing_samp))
    } else {
        worst_id <- NULL
    }
    
    # Remove worst reference if possible
    gene_data_subset <- gene_data_subset[, colnames(gene_data_subset) != worst_id]
    
}