#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# subset_core_gene.R smarcescens.tsv smarcescens_feat_seqs_renamed smarcescens.csv 30 300 smarcescens_core_genes smarcescens_feat_seqs
# args <- c('~/downloads/smarcescens.tsv',
#           '~/downloads/smarcescens_feat_seqs_renamed',
#           '~/downloads/smarcescens.csv',
#           '10', '300', '_no_group_defined__core_genes', '_no_group_defined__feat_seqs')
names(args) <- c("gene_families", "gene_seq_dir_path", "metadata", "min_core_genes",  "max_core_genes", "tsv_output_path", "fasta_output_path")
args <- as.list(args)
raw_gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = FALSE, col.names = c('sample_id', 'ref_id', 'ref_name', 'ref_desc', 'usage'), sep = '\t')
min_core_genes <- as.integer(args$min_core_genes)
max_core_genes <- as.integer(args$max_core_genes)

# Infer number of samples and references
total_count <- ncol(raw_gene_data) - 22
all_ids <- colnames(raw_gene_data)[23:ncol(raw_gene_data)]
sample_ids <- unique(metadata$sample_id[metadata$sample_id %in% all_ids])
ref_ids <- all_ids[! all_ids %in% sample_ids]

# Create matrix with number of genes for each gene in each sample
cluster_data <- raw_gene_data[, 23:ncol(raw_gene_data)]
cluster_data[, all_ids] <- lapply(cluster_data[, all_ids], function(column) {
    unlist(lapply(strsplit(column, split = '[;:]'), length))
})

# Create pairwise matrix of number of shared genes between all samples and references
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

# Make subset Pirate output for clusters
output_clusters <- lapply(clustered_samples, function(ids) {
    is_core_gene <- rowSums(present_and_single[, ids]) == length(ids)
    output <- raw_gene_data[is_core_gene, c(colnames(raw_gene_data)[1:22], ids)]
    if (nrow(output) > max_core_genes) {
        output <- output[1:max_core_genes, ]
    }
    return(output)
})
dir.create(args$tsv_output_path, showWarnings = FALSE)
for (index in seq_along(output_clusters)) {
    out_path <- file.path(args$tsv_output_path, paste0('cluster_', index, '.tsv'))
    write.table(output_clusters[[index]], file = out_path, row.names = FALSE, sep = '\t', quote = FALSE)
}

# Save the IDs of any samples and references that were removed from the analysis
removed_ids <- character()
# if (any(!has_ref_and_samp)) {
#     removed_ids <- c(removed_ids, unlist(clustered_samples[! has_ref_and_samp]))
# }
removed_sample_ids <- removed_ids[removed_ids %in% sample_ids]
removed_ref_ids <- removed_ids[removed_ids %in% ref_ids]
writeLines(removed_sample_ids, con = 'removed_sample_ids.txt')
writeLines(removed_ref_ids, con = 'removed_ref_ids.txt')

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

# Make directory to store subdirectories with each cluster's gene sequences
dir.create(args$fasta_output_path, showWarnings = FALSE)

# Save gene sequences for each cluster
for (index in seq_along(output_clusters)) {
    out_dir_path <- file.path(args$fasta_output_path, paste0('cluster_', index))
    dir.create(out_dir_path, showWarnings = FALSE)
    for (gene_id in output_clusters[[index]]$gene_family) {
        in_path <- file.path(args$gene_seq_dir_path, paste0(gene_id, '.fasta'))
        out_path <- file.path(out_dir_path, paste0(gene_id, '.fasta'))
        seqs <- read_fasta(in_path)
        names(seqs) <- trimws(names(seqs))
        samples_in_subset <- names(output_clusters[[index]])[names(output_clusters[[index]]) %in% c(sample_ids, ref_ids)]
        seqs <- seqs[samples_in_subset]
        seqs <- seqs[!is.na(seqs)]
        write_fasta(seqs, out_path)
    }
}


