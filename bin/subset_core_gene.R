#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- c('/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/2c/1b29f55ab50c5b8b16a6d19c9d2b52/subgroup.tsv',
          '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/2c/1b29f55ab50c5b8b16a6d19c9d2b52/subgroup_feat_seqs_renamed',
          '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/2c/1b29f55ab50c5b8b16a6d19c9d2b52/samplesheet.valid.csv',
          '10', '100', 'subgroup_core_genes', 'subgroup_feat_seqs')
names(args) <- c("gene_families", "gene_seq_dir_path", "metadata", "min_core_genes",  "max_core_genes", "csv_output_path", "fasta_output_path")
args <- as.list(args)
raw_gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = TRUE, sep = ',', row.names = NULL, check.names = FALSE)
min_core_genes <- as.integer(args$min_core_genes)
max_core_genes <- as.integer(args$max_core_genes)
raw_gene_data1 <- raw_gene_data

# Infer number of samples and references
total_count <- ncol(raw_gene_data) - 22
all_ids <- colnames(raw_gene_data)[23:ncol(raw_gene_data)]
sample_ids <- all_ids[all_ids %in% metadata$sample_id]
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Create matrix with number of genes for each gene in each sample
cluster_data <- raw_gene_data1[, 23:ncol(raw_gene_data1)]
cluster_data[, all_ids] <- lapply(cluster_data[, all_ids], function(column) {
  unlist(lapply(strsplit(column, split = '[;:]'), length))
})

# Create pairwise matrix of number of shared genes between all smaples and references
count_shared_genes <- function(data, sample_1, sample_2) {
    if (length(sample_1) != length(sample_2)) {
        stop("Arguments must be of same length.")
    }
    unlist(lapply(seq_along(sample_1), function(i) {
        if (sample_1[i] == sample_2[i]) {
            return(NA_integer_)
        } else {
            return(sum(data[[sample_1[i]]] == 1 & data[[sample_2[i]]] == 1))
        }
    }))
}
calc_shared_gene_matrix <- function(data, sample_ids) {
    names(sample_ids) <- sample_ids
    outer(sample_ids, sample_ids, function(a, b) count_shared_genes(data, a, b))
}
shared_genes <- calc_shared_gene_matrix(cluster_data, all_ids)

# Initialize "clusters" with one sample each
clustered_samples <- as.list(sample_ids)
names(clustered_samples) <- sample_ids # Sample IDs are also used for the initial cluster IDs

# Add the reference with the largest number of genes shared with each sample
#  This insures that at least 1 reference is included in each cluster
clustered_samples <- lapply(clustered_samples, function(id) {
    ref_id <- names(which.max(shared_genes[id, ]))
    n_shared <- shared_genes[id, ref_id]
    if (n_shared >= min_core_genes) {
        return(c(id, ref_id))
    } else {
        return(id)
    }
})

# Initialize shared gene matrix for clusters
present_and_single <- as.data.frame(lapply(cluster_data, function(x) x == 1), check.names = FALSE)

# Keep combining clusters until no more can be combined
# NOTE: there probably is a more efficient way of doing this.
#   This way was quick to code, but recalculates the shared genes between everything
#   each time clusters are combined even though most of the matrix does not change or is not used.
#   It would be faster to maintain and update a single table of shared gene counts for the clusters.
get_cluster_tf <- function(clusters) {
    as.data.frame(lapply(clusters, function(ids) {
        apply(present_and_single[, ids, drop = FALSE], 1, all)
    }), check.names = FALSE)
}
get_next_combination <- function(clusters, must_include = NULL) {
    cluster_tf <- get_cluster_tf(clusters)
    cluster_shared_genes <- calc_shared_gene_matrix(cluster_tf, names(clusters))
    if (!is.null(must_include)) {
        cluster_shared_genes <- cluster_shared_genes[must_include, , drop = FALSE]
    }
    max_shared <- max(cluster_shared_genes, na.rm = TRUE)
    max_shared_index <- which(cluster_shared_genes == max_shared, arr.ind = TRUE)
    to_combine <- c(rownames(cluster_shared_genes)[max_shared_index[1, "row"]],
                    colnames(cluster_shared_genes)[max_shared_index[1, "col"]])
    if (max_shared >= min_core_genes) {
        return(to_combine)
    } else {
        return(NULL)
    }
}
while (length(clustered_samples) > 1 && ! is.null(to_combine <- get_next_combination(clustered_samples))) {
    clustered_samples[[to_combine[1]]] <- unique(unlist(clustered_samples[to_combine]))
    clustered_samples[to_combine[2]] <- NULL
}

# Add as many references as possible to each cluster
# NOTE: there probably is a more efficient way of doing this.
#   This way was quick to code, but recalculates the shared genes between everything
#   each time clusters are combined even though most of the matrix does not change or is not used.
#   It would be faster to maintain and update a single table of shared gene counts for the clusters.
clustered_samples <- lapply(names(clustered_samples), function(cluster_id) {
    cluster_and_refs <- c(clustered_samples[cluster_id], setNames(ref_ids, ref_ids))
    while (length(cluster_and_refs) > 1 && ! is.null(to_combine <- get_next_combination(cluster_and_refs, must_include = cluster_id))) {
        cluster_and_refs[[to_combine[1]]] <- unique(unlist(cluster_and_refs[to_combine]))
        cluster_and_refs[to_combine[2]] <- NULL
    }
    return(cluster_and_refs[[1]])
})

# Determine which cluster have both samples and references
has_ref_and_samp <- unlist(lapply(clustered_samples, function(ids) {
    any(ids %in% ref_ids) && any(ids %in% sample_ids)
}))

# Make subset Pirate output for clusters with both samples and references
output_clusters <- lapply(clustered_samples[has_ref_and_samp], function(ids) {
    is_core_gene <- rowSums(present_and_single[, ids]) == length(ids)
    output <- raw_gene_data[is_core_gene, c(colnames(raw_gene_data)[1:22], ids)]
    if (nrow(output) > max_core_genes) {
        output <- output[1:max_core_genes, ]
    }
    return(output)
})
dir.create(args$csv_output_path, showWarnings = FALSE)
for (index in seq_along(output_clusters)) {
    out_path <- file.path(args$csv_output_path, paste0('cluster_', index, '.tsv'))
    write.table(output_clusters[[index]], file = out_path, row.names = FALSE, sep = '\t', quote = FALSE)
}

# Save the IDs of any samples and references that were removed from the analysis
if (any(!has_ref_and_samp)) {
    removed_ids <- unlist(clustered_samples[! has_ref_and_samp])
} else {
    removed_ids <- character()
}
removed_sample_ids <- removed_ids[removed_ids %in% sample_ids]
removed_ref_ids <- removed_ids[removed_ids %in% ref_ids]
writeLines(removed_sample_ids, con = 'removed_sample_ids.txt')
writeLines(removed_ref_ids, con = 'removed_ref_ids.txt')

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

# Make directory to store subdirectories with each cluster's gene sequences
dir.create(args$fasta_output_path, showWarnings = FALSE)

# Save gene sequences for each cluster
for (index in seq_along(output_clusters)) {
    out_dir_path <- file.path(args$fasta_output_path, paste0('cluster_', index))
    dir.create(out_dir_path, showWarnings = FALSE)
    for (gene_id in output_clusters[[index$gene_family) {
        in_path <- file.path(args$gene_seq_dir_path, paste0(gene_id, '.fasta'))
        out_path <- file.path(out_dir_path, paste0(gene_id, '.fasta'))
        seqs <- read_fasta(in_path)
        seqs <- seqs[c(sample_ids, ref_ids)]
        seqs <- seqs[!is.na(seqs)]
        write_fasta(seqs, out_path)
        paste(seqs, out_path)
    }
}


