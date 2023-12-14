#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c("PIRATE.gene_families.ordered.tsv", "metadata_medium.csv", "10", "0.8", "0.5") # NOTE: this is for testing, should be commented out for production code
args <- c("Bradyrhizobium_yuanmingense.tsv", "metadata_brady.csv", "10", "0.8", "0.5") # NOTE: this is for testing, should be commented out for production code
names(args) <- c("gene_families", "metadata", "min_core_genes", "min_core_samps", "min_core_refs")
args <- as.list(args)
raw_gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = TRUE, sep = ',')
min_core_genes <- as.integer(args$min_core_genes)

# Infer number of samples and references
total_count <- ncol(raw_gene_data) - 22
all_ids <- colnames(raw_gene_data)[23:ncol(raw_gene_data)]
sample_ids <- all_ids[all_ids %in% gsub(pattern = '[-.]+', replacement = '_', metadata$sample)] #TODO: change this once sample IDs are standardized
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Get minimum number of samples and references that is acceptable
min_core_samps <- ceiling(as.numeric(args$min_core_samps) * length(sample_ids))
min_core_refs <- ceiling(as.numeric(args$min_core_refs) * length(ref_ids))

# Replace gene name columns for each sample/ref with number of genes found
gene_data[, all_ids] <- lapply(raw_gene_data[, all_ids], function(column) {
    unlist(lapply(strsplit(column, split = '[;:]'), length))
})

# Remove samples until a core genome is found
gene_data_subset <- gene_data
current_sample_ids <- sample_ids
current_ref_ids <- ref_ids
get_n_core_single <- function(ids) {
    rowSums(gene_data_subset[, ids] == 1)
}
current_core_genes <- sum(get_n_core_single(current_sample_ids) == length(current_sample_ids) & 
                              get_n_core_single(current_ref_ids) == length(current_ref_ids))
while (current_core_genes < min_core_genes) {
    # Temporary debugging output TODO: delete when done
    print(paste0('samples: ', length(current_sample_ids)))
    print(paste0('refs:    ', length(current_ref_ids)))
    print(paste0('core:    ', current_core_genes))
    
    # Find how many genes each sample are missing or multi-copy
    n_bad_samp <- colSums(gene_data_subset[, current_sample_ids] != 1)
    n_bad_ref <- colSums(gene_data_subset[, current_ref_ids] != 1)
    
    # Remove worst sample/ref if possible, preferring to remove references
    if (length(current_sample_ids) > min_core_samps && length(current_ref_ids) > min_core_refs) {
        worst_id <- names(which.min(c(n_bad_samp, n_bad_ref)))
    } else if (length(current_ref_ids) > min_core_refs) {
        worst_id <- names(which.min(n_bad_ref))
    } else if (length(current_sample_ids) > min_core_samps) {
        worst_id <- names(which.min(n_bad_samp))
    } else {
        worst_id <- NULL
        break
    }
    
    # Remove worst sample/reference
    gene_data_subset <- gene_data_subset[, colnames(gene_data_subset) != worst_id]
    current_sample_ids <- current_sample_ids[current_sample_ids != worst_id]
    current_ref_ids <- current_ref_ids[current_ref_ids != worst_id]

    current_core_genes <- sum(get_n_core_single(current_sample_ids) == length(current_sample_ids) & 
                                  get_n_core_single(current_ref_ids) == length(current_ref_ids))
}


# Filter out non-core genes 
output <- raw_gene_data[get_n_core_single(current_sample_ids) == length(current_sample_ids) & 
                            get_n_core_single(current_ref_ids) == length(current_ref_ids), ]

# Write filtered output table
write.table(output, file = 'test_output.tsv', row.names = FALSE, sep = '\t', quote = FALSE)
