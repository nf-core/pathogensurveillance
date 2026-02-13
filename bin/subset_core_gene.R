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
message_data_path <- 'message_data.tsv'
tsv_output_path <- 'core_genes'
fasta_output_path <- 'feat_seqs'

# Initialize message data to store messages shown to the user
message_data <- data.frame(
    sample_id = character(0),
    report_group_id = character(0),
    reference_id = character(0),
    workflow = character(0),
    message_type = character(0),
    description = character(0)
)

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
names(args) <- c("gene_families", "gene_seq_dir_path", "metadata", "min_genes",  "max_genes", "prefix")
args <- as.list(args)
raw_gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = FALSE, col.names = c('sample_id', 'ref_id', 'ref_name', 'ref_desc', 'usage'), sep = '\t')
min_genes <- as.integer(args$min_genes)
max_genes <- as.integer(args$max_genes)

# Replace runs of underscores with a single underscore in sample IDs since this is what PIRATE does
metadata$modified_sample_ids <- gsub(metadata$sample_id, pattern = '_+', replacement = '_')

# Infer number of samples and references
total_count <- ncol(raw_gene_data) - 22
all_ids <- colnames(raw_gene_data)[23:ncol(raw_gene_data)]
sample_ids <- unique(metadata$modified_sample_ids[metadata$modified_sample_ids %in% all_ids])
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Create matrix with TRUE/FALSE for whether a gene is present and single copy
present_and_single_original <- as.data.frame(lapply(raw_gene_data[, all_ids], function(column) {
    unlist(lapply(strsplit(column, split = '[;:]'), length)) == 1
}), check.names = FALSE)
present_and_single <- present_and_single_original

# Initialize "clusters" with one sample each
# Sample IDs are also used for the cluster IDs, although the specific sample used is arbitrary
clusters <- as.list(all_ids)
names(clusters) <- all_ids

# Make table with number of genes shared for all pairwise comparisons
count_shared_genes <- function(cluster_1_id, cluster_2_id) {
    count <- sum(present_and_single[[cluster_1_id]] & present_and_single[[cluster_2_id]])
    data.frame(cluster_1_id = cluster_1_id, cluster_2_id = cluster_2_id, count = count)
}
cluster_shared_count <- do.call(rbind, lapply(combn(seq_along(all_ids), 2, simplify = F), function(x) count_shared_genes(all_ids[x[1]], all_ids[x[2]])))

# Keep combining clusters until no more can be combined
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
while (length(clusters) > 1 && ! is.null(best_pair <- get_next_combination())) {
    # Update sample clusters
    best_ids <- c(best_pair$cluster_1_id, best_pair$cluster_2_id)
    clusters[[best_pair$cluster_1_id]] <- unname(unique(unlist(clusters[best_ids])))
    clusters[best_pair$cluster_2_id] <- NULL

    # Update gene presence/absence matrix
    present_and_single[[best_pair$cluster_1_id]] <- present_and_single[[best_pair$cluster_1_id]] & present_and_single[[best_pair$cluster_2_id]]
    present_and_single[best_pair$cluster_2_id] <- NULL

    # Update pairwise shared gene counts between clusters
    cluster_shared_count <- cluster_shared_count[! (cluster_shared_count$cluster_1_id %in% best_ids | cluster_shared_count$cluster_2_id %in% best_ids), ]
    other_clusters <- names(clusters)[names(clusters) != best_pair$cluster_1_id]
    new_shared_count <- do.call(rbind, lapply(other_clusters, function(n) count_shared_genes(n, best_pair$cluster_1_id)))
    cluster_shared_count <- rbind(cluster_shared_count, new_shared_count)

    # Save updated clusters and cluster statistics
    all_clusterings <- c(all_clusterings, list(clusters))
    clustering_stats <- c(clustering_stats, list(calc_cluster_stats(clusters)))
}

# Choose best clustering based on clustering statistics
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
clustering_stats$mean_cluster_size <- ifelse(clustering_stats$valid_cluster_count == 0, 0, clustering_stats$sample_count / clustering_stats$valid_cluster_count)
if (max(clustering_stats$mean_cluster_size, na.rm = T) == 0) {
  clustering_stats$mean_cluster_size_score <- 0
} else {
  clustering_stats$mean_cluster_size_score <- clustering_stats$mean_cluster_size / max(clustering_stats$mean_cluster_size, na.rm = T)
}
non_zero_cluster_range <-  range(clustering_stats$valid_cluster_count[clustering_stats$valid_cluster_count != 0])
clustering_stats$cluster_count_score <- ifelse(clustering_stats$valid_cluster_count == 0, 0, 1 / (clustering_stats$valid_cluster_count - min(non_zero_cluster_range) + 1))
clustering_stats$shared_gene_score <- ifelse(clustering_stats$weigthed_shared_genes >= max_genes, 1, logistic_scaling_func(clustering_stats$weigthed_shared_genes))
clustering_stats$overall_score <- vapply(seq_len(nrow(clustering_stats)), FUN.VALUE = numeric(1), function(i) {
    sum(vapply(names(clustering_weights), FUN.VALUE = numeric(1), function(n) clustering_stats[[n]][i] * clustering_weights[n]))
})
best_clusters <- all_clusterings[[which.max(clustering_stats$overall_score)]]

# Remove clusters that have too few samples/references
is_valid <- vapply(best_clusters, is_valid_cluster, FUN.VALUE = logical(1))
removed_ids <- unlist(best_clusters[! is_valid])
removed_sample_ids <- removed_ids[removed_ids %in% sample_ids]
removed_reference_ids <- removed_ids[removed_ids %in% ref_ids]
best_clusters <- best_clusters[is_valid]

# Report removed samples and references
report_group_id <- gsub(basename(args$metadata), pattern = '\\.tsv$', replacement = '')
if (length(removed_ids) > 0) {
    warning('Removed ', sum(! is_valid), ' clusters with fewer than ', min_cluster_size, ' samples, totaling ',
            length(removed_sample_ids), ' samples and ', length(removed_reference_ids), ' references.')
}
if (length(removed_sample_ids) > 0) {
    message_data <- rbind(message_data, data.frame(
        sample_id = removed_sample_ids,
        report_group_id = report_group_id,
        reference_id = NA_character_,
        workflow = 'CORE_GENOME_PHYLOGENY',
        message_type = 'WARNING',
        description = 'Sample removed from core gene phylogeny in order to find enough core genes among other samples/references.'
    ))
}
if (length(removed_reference_ids) > 0) {
    message_data <- rbind(message_data, data.frame(
        sample_id = NA_character_,
        report_group_id = report_group_id,
        reference_id = removed_reference_ids,
        workflow = 'CORE_GENOME_PHYLOGENY',
        message_type = 'NOTE',
        description = 'Reference removed from core gene phylogeny in order to find enough core genes among other samples/references.'
    ))
}

# Make subset Pirate output for clusters
output_cluster_data <- lapply(best_clusters, function(ids) {
    is_core_gene <- rowSums(present_and_single_original[, ids]) == length(ids)
    output <- raw_gene_data[is_core_gene, c(colnames(raw_gene_data)[1:22], ids)]
    if (nrow(output) > max_genes) {
        output <- output[1:max_genes, ]
    }
    return(output)
})
dir.create(tsv_output_path, showWarnings = FALSE)
for (index in seq_along(output_cluster_data)) {
    out_path <- file.path(tsv_output_path, paste0(args$prefix, '--cluster_', index, '.tsv'))
    write.table(output_cluster_data[[index]], file = out_path, row.names = FALSE, sep = '\t', quote = FALSE)
}

# Copy sequences for selected genes and samples and put in new folder
read_fasta <- function(file_path) {
    # Read raw string
    raw_data <- paste0(readLines(file_path), collapse = '\n')

    # Return an empty vector an a warning if no sequences are found
    if (raw_data == "") {
        warning(paste0("No sequences found in the file: ", file_path))
        return(character(0))
    }

    # Find location of every header start
    split_data <- strsplit(raw_data, "\n>")[[1]]

    # Split the data for each sequence into lines
    split_data <- strsplit(split_data, split = "\n")

    # The first lines are headers, so remove those
    headers <- vapply(split_data, FUN = `[`, FUN.VALUE = character(1), 1)
    split_data <- lapply(split_data, FUN = `[`, -1)

    # Remove the > from the first sequence. The others were removed by the split
    headers[1] <- sub(headers[1], pattern = "^>", replacement = "")

    # Combine multiple lines into single sequences
    seqs <- vapply(split_data, FUN = paste0, FUN.VALUE = character(1), collapse = "")

    # Combine and return results
    return(stats::setNames(seqs, headers))
}
write_fasta <- function(seqs, path) {
    output <- paste0('>', names(seqs), '\n', seqs)
    writeLines(output, path)
}
dir.create(fasta_output_path, showWarnings = FALSE)
for (index in seq_along(output_cluster_data)) {
    out_dir_path <- file.path(fasta_output_path, paste0(args$prefix, '--cluster_', index))
    dir.create(out_dir_path, showWarnings = FALSE)
    for (gene_id in output_cluster_data[[index]]$gene_family) {
        in_path <- file.path(args$gene_seq_dir_path, paste0(gene_id, '.fasta'))
        out_path <- file.path(out_dir_path, paste0(gene_id, '.fasta'))
        seqs <- read_fasta(in_path)
        names(seqs) <- trimws(names(seqs))
        names(seqs) <- sub(names(seqs), pattern = '.*genome:(.*)\tgene:.*', replacement = '\\1')
        samples_in_subset <- names(output_cluster_data[[index]])[names(output_cluster_data[[index]]) %in% c(sample_ids, ref_ids)]
        seqs <- seqs[samples_in_subset]
        seqs <- seqs[!is.na(seqs)]
        write_fasta(seqs, out_path)
    }
}

# Write data for messages to be shown to the user, such as warnings about removed samples
write.table(message_data, file = message_data_path, row.names = FALSE, na = '', sep = '\t')

# Save clustering statistics
write.table(clustering_stats, file = clustering_stats_path, row.names = FALSE, na = '', sep = '\t')
