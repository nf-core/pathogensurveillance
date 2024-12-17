#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     "~/projects/pathogensurveillance/work/9c/8c88175fd78b3128ef41661497227a/_no_group_defined_.csv",
#     "10",
#     "300",
#     "deleteme",
#     "~/projects/pathogensurveillance/work/9c/8c88175fd78b3128ef41661497227a/GCF_027943255_1-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/9c/8c88175fd78b3128ef41661497227a/GCA_005966545_1-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/9c/8c88175fd78b3128ef41661497227a/GCA_026225685_1-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/9c/8c88175fd78b3128ef41661497227a/SRR29399310-eukaryota_odb10-busco",
#     "~/projects/pathogensurveillance/work/9c/8c88175fd78b3128ef41661497227a/GCA_026401395_1-eukaryota_odb10-busco"
# )

args <- as.list(args)
if (length(args) < 5) {
    stop('No busco result directories supplied.')
}
busco_dirs <- unlist(args[5:length(args)])
args <- args[1:4]
names(args) <- c("sample_ref_pairs", "min_genes", "max_genes",  "feat_seqs_out_path")
metadata <- read.csv(args$sample_ref_pairs, header = FALSE, col.names = c('sample_id', 'ref_id', 'ref_name', 'ref_desc', 'usage'))
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

# Create pairwise matrix of number of shared genes between all samples and references
count_shared_genes <- function(data, sample_1, sample_2) {
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
    if (n_shared >= min_genes) {
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
    if (max_shared >= min_genes) {
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
clustered_samples <- clustered_samples[has_ref_and_samp]

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

# Make directory to store subdirectories with each cluster's gene sequences
dir.create(args$feat_seqs_out_path, showWarnings = FALSE)

# Save gene sequences for each cluster
output_clusters <- lapply(clustered_samples, function(ids) {
    is_core_gene <- rowSums(present_and_single[, ids]) == length(ids)
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


