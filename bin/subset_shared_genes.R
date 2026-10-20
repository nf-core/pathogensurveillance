#!/usr/bin/env Rscript

# MIT License
#
# Copyright (c) Zachary S.L. Foster and Niklaus J. Grunwald
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Subset shared genes from a unified gene_sequences/ directory and prepare
# per-cluster FASTA files for alignment and tree building.
# Args: <gene_seq_dir> <samp_ref_pairs> <min_genes> <max_genes> <prefix>

# Constants
min_cluster_size <- 4
clustering_weights <- c(
    sample_proportion = 5,
    reference_proportion = 2,
    shared_gene_score = 1,
    mean_cluster_size_score = 1
)
clustering_weights <- clustering_weights / sum(clustering_weights)
clustering_stats_path <- 'clustering_statistics.tsv'

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c("gene_seq_dir", "samp_ref_pairs", "min_genes", "max_genes", "prefix")
args <- as.list(args)
min_genes <- as.integer(args$min_genes)
max_genes <- as.integer(args$max_genes)
prefix <- args$prefix

metadata <- read.csv(args$samp_ref_pairs, header = FALSE,
                     col.names = c('sample_id', 'ref_id', 'ref_name', 'ref_desc', 'usage'),
                     sep = '\t', stringsAsFactors = FALSE)
# PIRATE collapses runs of underscores in genome IDs, so match using the same IDs
metadata$sample_id <- gsub('_+', '_', metadata$sample_id)

# ---------------------------------------------------------------------------
# Read gene sequences and build presence/absence matrix
# ---------------------------------------------------------------------------
read_fasta <- function(file_path) {
    raw_data <- paste0(readLines(file_path), collapse = '\n')
    if (raw_data == "") {
        return(character(0))
    }
    split_data <- strsplit(raw_data, "\n>")[[1]]
    split_data <- strsplit(split_data, split = "\n")
    headers <- vapply(split_data, FUN = `[`, FUN.VALUE = character(1), 1)
    split_data <- lapply(split_data, FUN = `[`, -1)
    headers[1] <- sub(headers[1], pattern = "^>", replacement = "")
    seqs <- vapply(split_data, FUN = paste0, FUN.VALUE = character(1), collapse = "")
    names(seqs) <- trimws(headers)
    return(seqs)
}

gene_files <- list.files(args$gene_seq_dir, pattern = "\\.fasta$", full.names = TRUE)
if (length(gene_files) == 0) {
    stop('No gene fasta files found in ', args$gene_seq_dir)
}

gene_ids <- sub("\\.fasta$", "", basename(gene_files))
names(gene_files) <- gene_ids

# Collect all genome IDs
all_ids <- unique(unlist(lapply(gene_files, function(f) {
    names(read_fasta(f))
})))
if (length(all_ids) == 0) {
    stop('No sequences found in gene files.')
}

sample_ids <- unique(metadata$sample_id[metadata$sample_id %in% all_ids])
ref_ids <- all_ids[! all_ids %in% sample_ids]

if (length(sample_ids) == 0) {
    stop('No samples are present in gene sequence data.')
}

# Build presence/absence matrix: TRUE iff genome appears exactly once for a gene
# Rows are genes and columns are genome IDs so clusters can be manipulated as columns
present_matrix <- matrix(FALSE, nrow = length(gene_ids), ncol = length(all_ids),
                         dimnames = list(gene_ids, all_ids))
for (i in seq_along(gene_files)) {
    seqs <- read_fasta(gene_files[i])
    genome_counts <- table(names(seqs))
    present_matrix[i, ] <- vapply(all_ids, FUN.VALUE = logical(1), function(g) {
        !is.na(genome_counts[g]) && genome_counts[g] == 1
    })
}
present_and_single_original <- as.data.frame(present_matrix, check.names = FALSE)
present_and_single <- present_and_single_original

# ---------------------------------------------------------------------------
# Hierarchical clustering of genomes by shared single-copy genes
# ---------------------------------------------------------------------------
clusters <- as.list(all_ids)
names(clusters) <- all_ids

count_shared_genes <- function(cluster_1_id, cluster_2_id) {
    count <- sum(present_and_single[[cluster_1_id]] & present_and_single[[cluster_2_id]])
    data.frame(cluster_1_id = cluster_1_id, cluster_2_id = cluster_2_id, count = count)
}

cluster_shared_count <- do.call(rbind, lapply(combn(seq_along(all_ids), 2, simplify = FALSE),
                                              function(x) count_shared_genes(all_ids[x[1]], all_ids[x[2]])))

get_next_combination <- function() {
    most_shared <- which.max(cluster_shared_count$count)
    if (cluster_shared_count$count[most_shared] >= min_genes) {
        return(cluster_shared_count[most_shared, ])
    } else {
        return(NULL)
    }
}

is_valid_cluster <- function(cluster) {
    length(cluster) >= min_cluster_size && any(cluster %in% sample_ids)
}

calc_cluster_stats <- function(clusters) {
    cluster_lengths <- vapply(clusters, length, FUN.VALUE = numeric(1))
    shared_gene_counts <- vapply(names(clusters), FUN.VALUE = numeric(1), function(id) {
        sum(present_and_single[[id]])
    })
    is_valid <- vapply(clusters, FUN.VALUE = logical(1), is_valid_cluster)
    data.frame(
        cluster_count = length(clusters),
        valid_cluster_count = sum(is_valid),
        sample_count = sum(unlist(clusters[is_valid]) %in% sample_ids),
        ref_count = sum(unlist(clusters[is_valid]) %in% ref_ids),
        weigthed_shared_genes = sum(shared_gene_counts * cluster_lengths / sum(cluster_lengths))
    )
}

all_clusterings <- list(clusters)
clustering_stats <- list(calc_cluster_stats(clusters))
while (length(clusters) > 1 && !is.null(best_pair <- get_next_combination())) {
    best_ids <- c(best_pair$cluster_1_id, best_pair$cluster_2_id)
    clusters[[best_pair$cluster_1_id]] <- unname(unique(unlist(clusters[best_ids])))
    clusters[best_pair$cluster_2_id] <- NULL

    present_and_single[[best_pair$cluster_1_id]] <- present_and_single[[best_pair$cluster_1_id]] & present_and_single[[best_pair$cluster_2_id]]
    present_and_single[best_pair$cluster_2_id] <- NULL

    cluster_shared_count <- cluster_shared_count[!(cluster_shared_count$cluster_1_id %in% best_ids | cluster_shared_count$cluster_2_id %in% best_ids), ]
    other_clusters <- names(clusters)[names(clusters) != best_pair$cluster_1_id]
    new_shared_count <- do.call(rbind, lapply(other_clusters, function(n) count_shared_genes(n, best_pair$cluster_1_id)))
    cluster_shared_count <- rbind(cluster_shared_count, new_shared_count)

    all_clusterings <- c(all_clusterings, list(clusters))
    clustering_stats <- c(clustering_stats, list(calc_cluster_stats(clusters)))
}

# ---------------------------------------------------------------------------
# Choose best clustering
# ---------------------------------------------------------------------------
logistic_scaling_func <- function(shared_count, floor = 0.1, ceiling = 1, midpoint = 0, steepness = 0.02) {
    logistic_value <- ceiling / (1 + exp(1)^(-steepness * (shared_count - midpoint)))
    (logistic_value - 0.5 + floor) / (0.5 + floor)
}

clustering_stats <- do.call(rbind, clustering_stats)
clustering_stats$sample_proportion <- clustering_stats$sample_count / length(sample_ids)
if (length(ref_ids) == 0) {
    clustering_stats$reference_proportion <- 1
} else {
    clustering_stats$reference_proportion <- clustering_stats$ref_count / length(ref_ids)
}
clustering_stats$mean_cluster_size <- ifelse(clustering_stats$valid_cluster_count == 0, 0,
                                             clustering_stats$sample_count / clustering_stats$valid_cluster_count)
if (max(clustering_stats$mean_cluster_size, na.rm = TRUE) == 0) {
    clustering_stats$mean_cluster_size_score <- 0
} else {
    clustering_stats$mean_cluster_size_score <- clustering_stats$mean_cluster_size / max(clustering_stats$mean_cluster_size, na.rm = TRUE)
}
non_zero_cluster_range <- range(clustering_stats$valid_cluster_count[clustering_stats$valid_cluster_count != 0])
clustering_stats$cluster_count_score <- ifelse(clustering_stats$valid_cluster_count == 0, 0,
                                               1 / (clustering_stats$valid_cluster_count - min(non_zero_cluster_range) + 1))
clustering_stats$shared_gene_score <- ifelse(clustering_stats$weigthed_shared_genes >= max_genes, 1,
                                             logistic_scaling_func(clustering_stats$weigthed_shared_genes))
clustering_stats$overall_score <- vapply(seq_len(nrow(clustering_stats)), FUN.VALUE = numeric(1), function(i) {
    sum(vapply(names(clustering_weights), FUN.VALUE = numeric(1), function(n) clustering_stats[[n]][i] * clustering_weights[n]))
})
best_clusters <- all_clusterings[[which.max(clustering_stats$overall_score)]]

# Remove invalid clusters and record removed IDs
is_valid <- vapply(best_clusters, is_valid_cluster, FUN.VALUE = logical(1))
removed_ids <- unlist(best_clusters[!is_valid])
removed_sample_ids <- removed_ids[removed_ids %in% sample_ids]
removed_reference_ids <- removed_ids[removed_ids %in% ref_ids]
best_clusters <- best_clusters[is_valid]

removed_id_data <- data.frame(
    id = c(removed_sample_ids, removed_reference_ids),
    type = c(rep('sample', length(removed_sample_ids)), rep('reference', length(removed_reference_ids))),
    stringsAsFactors = FALSE
)
write.table(removed_id_data, file = 'removed_ids.tsv', sep = '\t', quote = FALSE, row.names = FALSE)

# ---------------------------------------------------------------------------
# Write per-cluster gene sequences and gene lists
# ---------------------------------------------------------------------------
feat_seqs_dir <- 'feat_seqs'
dir.create(feat_seqs_dir, showWarnings = FALSE)

for (index in seq_along(best_clusters)) {
    cluster_members <- best_clusters[[index]]
    is_core_gene <- rowSums(present_and_single_original[, cluster_members, drop = FALSE]) == length(cluster_members)
    cluster_gene_ids <- gene_ids[is_core_gene]
    if (length(cluster_gene_ids) > max_genes) {
        cluster_gene_ids <- cluster_gene_ids[1:max_genes]
    }

    cluster_name <- paste0(prefix, '--cluster_', index)
    cluster_dir <- file.path(feat_seqs_dir, cluster_name)
    dir.create(cluster_dir, recursive = TRUE, showWarnings = FALSE)

    for (gene_id in cluster_gene_ids) {
        seqs <- read_fasta(gene_files[gene_id])
        seqs <- seqs[names(seqs) %in% cluster_members]
        seqs <- seqs[!is.na(seqs)]
        if (length(seqs) == 0) next

        out_path <- file.path(cluster_dir, paste0(cluster_name, '--', gene_id, '.fasta'))
        writeLines(paste0('>', names(seqs), '\n', seqs), out_path)
    }

    # Write gene list for this cluster
    gene_list_path <- paste0(cluster_name, '_gene_list.tsv')
    write.table(data.frame(gene_name = cluster_gene_ids),
                file = gene_list_path, sep = '\t', quote = FALSE, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# Write clustering statistics
# ---------------------------------------------------------------------------
write.table(clustering_stats, file = clustering_stats_path, row.names = FALSE, na = '', sep = '\t')

cat("subset_shared_genes:\n  r-base: ", R.version.string, "\n", sep = "")
