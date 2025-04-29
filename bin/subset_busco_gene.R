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
min_cluster_size <- 3
clustering_weights <- c(
    sample_proportion = 5,
    reference_proportion = 2,
    shared_gene_score = 1,
    mean_cluster_size_score = 1
)
clustering_weights <- clustering_weights / sum(clustering_weights)
clustering_stats_path <- 'clustering_statistics.tsv'
message_data_path <- 'message_data.tsv'

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
# args <- c(
#     "~/projects/pathogensurveillance/work/3c/39705646a22b25e90ae0e88a7bd931/_no_group_defined_.tsv",
#     "10",
#     "300",
#     "deleteme",
#     "~/projects/pathogensurveillance/work/3c/39705646a22b25e90ae0e88a7bd931/GCF_000002765_6-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/3c/39705646a22b25e90ae0e88a7bd931/SRR27942447-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/3c/39705646a22b25e90ae0e88a7bd931/SRR29397006-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/3c/39705646a22b25e90ae0e88a7bd931/SRR29695256-eukaryota_odb10-busco"
# )
args <- as.list(args)
if (length(args) < 5) {
    stop('No busco result directories supplied.')
}
busco_dirs <- unlist(args[5:length(args)])
args <- args[1:4]
names(args) <- c("metadata", "min_genes", "max_genes",  "feat_seqs_out_path")
metadata <- read.csv(args$metadata, header = FALSE, col.names = c('sample_id', 'ref_id', 'ref_name', 'ref_desc', 'usage'), sep = '\t')
min_genes <- as.integer(args$min_genes)
max_genes <- as.integer(args$max_genes)

# Figure out which ref/sample id goes with which busco directory
file_path_wth_db_id <- list.files(busco_dirs[1], recursive = TRUE, include.dirs = TRUE, pattern = '^run_')
lineage_db_id <- sub(basename(file_path_wth_db_id), pattern = '^run_', replacement = '')
names(busco_dirs) <- sub(basename(busco_dirs), pattern = paste0('-', lineage_db_id, '-busco$'), replacement = '')

# Infer number of samples and references
all_ids <- names(busco_dirs)
sample_ids <- unique(metadata$sample_id[metadata$sample_id %in% all_ids])
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Create matrix with number of genes for each gene in each sample
file_name_without_ext <- function(paths, ext) {
    sub(basename(paths), pattern = paste0('\\.', ext,'$'), replacement = '')
}
gene_seq_paths <- lapply(busco_dirs, function(dir_path) {
    list.files(file.path(dir_path), recursive = TRUE, pattern = '.fna')
})
gene_seq_ids <- unique(file_name_without_ext(unlist(gene_seq_paths), 'fna'))
cluster_data <- do.call(rbind, lapply(gene_seq_paths, function(paths) {
    single_copy_ids <- file_name_without_ext(paths[grepl(paths, pattern = 'single_copy_busco_sequences')], 'fna')
    multi_copy_ids <- file_name_without_ext(paths[grepl(paths, pattern = 'multi_copy_busco_sequences')], 'fna')
    fragmented_ids <- file_name_without_ext(paths[grepl(paths, pattern = 'fragmented_busco_sequences')], 'fna')
    out <- as.data.frame(as.list(stats::setNames(rep(0, length(gene_seq_ids)), gene_seq_ids)), check.names = FALSE)
    out[1, single_copy_ids] <- 1
    out[1, multi_copy_ids] <- 0
    out[1, fragmented_ids] <- 0
    return(out)
}))
cluster_data <- as.data.frame(t(cluster_data))

# Create matrix with TRUE/FALSE for whether a gene is present and single copy
present_and_single_original <- as.data.frame(lapply(cluster_data, function(column) {
    column == 1
}))
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
clustering_stats$mean_cluster_size_score <- clustering_stats$mean_cluster_size / max(clustering_stats$mean_cluster_size, na.rm = T)
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
        description = 'Sample removed from core gene phylogeny in order to find enough core genes amoung other samples/references.'
    ))
}
if (length(removed_reference_ids) > 0) {
    message_data <- rbind(message_data, data.frame(
        sample_id = NA_character_,
        report_group_id = report_group_id,
        reference_id = removed_reference_ids,
        workflow = 'CORE_GENOME_PHYLOGENY',
        message_type = 'NOTE',
        description = 'Reference removed from core gene phylogeny in order to find enough core genes amoung other samples/references.'
    ))
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

    # The first lines are headers, so remvove those
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
dir.create(args$feat_seqs_out_path, showWarnings = FALSE)
output_clusters <- lapply(best_clusters, function(ids) {
    is_core_gene <- rowSums(present_and_single_original[, ids]) == length(ids)
    output <- cluster_data[is_core_gene, ids]
    if (nrow(output) > max_genes) {
        output <- output[1:max_genes, ]
    }
    return(output)
})
for (index in seq_along(output_clusters)) {
    out_dir_path <- file.path(args$feat_seqs_out_path, paste0('cluster_', index))
    gene_ids <- rownames(output_clusters[[index]])
    sample_ids <- colnames(output_clusters[[index]])
    dir.create(out_dir_path, showWarnings = FALSE)
    for (gene_id in gene_ids) {
        seq_paths <- vapply(busco_dirs[sample_ids], list.files, pattern = paste0(gene_id, '\\.fna'), recursive = TRUE, full.names = TRUE, FUN.VALUE = character(1))
        seqs <- unlist(lapply(seq_along(sample_ids), function(i) {
            out <- read_fasta(seq_paths[i])
            names(out) <- sample_ids[i]
            return(out)
        }))
        out_path <- file.path(out_dir_path, paste0(gene_id, '.fasta'))
        write_fasta(seqs, out_path)
    }
}

# Write data for messages to be shown to the user, such as warnings about removed samples
write.table(message_data, file = message_data_path, row.names = FALSE, na = '', sep = '\t')

# Save clustering statistics
write.table(clustering_stats, file = clustering_stats_path, row.names = FALSE, na = '', sep = '\t')


