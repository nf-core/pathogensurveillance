#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c('/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/05/742d105a2c32516298767879e33bb4/subgroup.tsv',
#           '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/05/742d105a2c32516298767879e33bb4/subgroup_feat_seqs_renamed',
#           '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/05/742d105a2c32516298767879e33bb4/samplesheet.valid.csv',
#           '10', '0.8', '0.5', '100', 'subgroup_core_genes.tsv', 'subgroup_feat_seqs')
names(args) <- c("gene_families", "gene_seq_dir_path", "metadata", "min_core_genes", "min_core_samps", "min_core_refs", "max_core_genes", "csv_output_path", "fasta_output_path")
args <- as.list(args)
raw_gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = TRUE, sep = ',', row.names = NULL, check.names = FALSE)
min_core_genes <- as.integer(args$min_core_genes)
max_core_genes <- as.integer(args$max_core_genes)
min_core_samps <- as.integer(args$min_core_samps)
min_core_refs <- as.integer(args$min_core_refs)

# Infer number of samples and references
total_count <- ncol(raw_gene_data) - 22
all_ids <- colnames(raw_gene_data)[23:ncol(raw_gene_data)]
sample_ids <- all_ids[all_ids %in% metadata$sample_id]
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Get minimum number of samples and references that is acceptable
min_core_samps <- ceiling(min_core_samps * length(sample_ids))
min_core_refs <- ceiling(min_core_refs * length(ref_ids))

# Remove rows that cannot meet the minimum number of genomes
raw_gene_data <- raw_gene_data[raw_gene_data$number_genomes >= min_core_samps + min_core_refs, ]

# Replace gene name columns for each sample/ref with number of genes found
gene_data <- raw_gene_data
gene_data[, all_ids] <- lapply(raw_gene_data[, all_ids], function(column) {
    unlist(lapply(strsplit(column, split = '[;:]'), length))
})

# Remove samples until a core genome is found
gene_data_subset <- gene_data
current_sample_ids <- sample_ids
current_ref_ids <- ref_ids
get_n_core_single <- function(ids) {
    rowSums(gene_data_subset[, ids, drop = FALSE] == 1)
}
current_core_genes <- sum(get_n_core_single(current_sample_ids) == length(current_sample_ids) &
                              get_n_core_single(current_ref_ids) == length(current_ref_ids))
while (current_core_genes < min_core_genes) {
    # Temporary debugging output TODO: delete when done
    # print(paste0('samples:  ', length(current_sample_ids)))
    # print(paste0('refs:     ', length(current_ref_ids)))
    # print(paste0('core:     ', current_core_genes))

    # Find how many genes each sample are missing or multi-copy
    n_bad_samp <- colSums(gene_data_subset[, current_sample_ids] != 1)
    n_bad_ref <- colSums(gene_data_subset[, current_ref_ids] != 1)

    # Remove worst sample/ref if possible, preferring to remove references
    if (length(current_sample_ids) > min_core_samps && length(current_ref_ids) > min_core_refs) {
        worst_id <- names(which.max(c(n_bad_samp, n_bad_ref)))
    } else if (length(current_ref_ids) > min_core_refs) {
        worst_id <- names(which.max(n_bad_ref))
    } else if (length(current_sample_ids) > min_core_samps) {
        worst_id <- names(which.max(n_bad_samp))
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

# Remove excess genes if more than needed were found
if (nrow(output) > max_core_genes) {
    output <- output[1:max_core_genes, ]
}

# Write filtered output table
write.table(output, file = args$csv_output_path, row.names = FALSE, sep = '\t', quote = FALSE)

# Copy sequences for selected genes and samples and put in new folder
read_fasta <- function(path) {
    # Read file as a single character
    raw <- paste0(base::readLines(path), collapse = "\n")

    # Split by header
    raw_seqs <- strsplit(raw, split = '\n>')[[1]]

    # Split header and sequence
    split_seqs <- strsplit(raw_seqs, split = '\n')

    # Format as character vector named by header
    seqs <- unlist(lapply(split_seqs, function(x) gsub(x[2], pattern = '\n', replacement = '')))
    names(seqs) <- unlist(lapply(split_seqs, function(x) { # Some headers have tabs separated values and some just have the ID. Not sure if that's a bug in a previous step, but it is handled here anyway
        trimws(strsplit(split = '\t', x[1])[[1]][1])
    }))
    names(seqs)[1] <- sub(names(seqs)[1], pattern = '^>', replacement = '')

    return(seqs)
}

write_fasta <- function(seqs, path) {
    output <- paste0('>', names(seqs), '\n', seqs)
    writeLines(output, path)
}

dir.create(args$fasta_output_path, showWarnings = FALSE)
passing_sample_ids <- c(current_sample_ids, current_ref_ids)
for (gene_id in output$gene_family) {
    in_path <- file.path(args$gene_seq_dir_path, paste0(gene_id, '.fasta'))
    out_path <- file.path(args$fasta_output_path, paste0(gene_id, '.fasta'))
    seqs <- read_fasta(in_path)
    seqs <- seqs[passing_sample_ids]
    seqs <- seqs[!is.na(seqs)]
    write_fasta(seqs, out_path)
}

