#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- c('subset_core_inputs_for_logan/subgroup.tsv',
          'subset_core_inputs_for_logan/subgroup_feat_seqs_renamed',
          'subset_core_inputs_for_logan/samplesheet.valid.csv',
          '10', '100', 'deleteme/subgroup_core_genes', 'deleteme/subgroup_feat_seqs')
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

#Splitting trees if specified
cluster_data <- raw_gene_data1[, 23:ncol(raw_gene_data1)]
cluster_data[, all_ids] <- lapply(cluster_data[, all_ids], function(column) {
  unlist(lapply(strsplit(column, split = '[;:]'), length))
})

# treating multiple copies as missing
cluster_data[cluster_data > 1] <- 0

# Samples that don't have enough genes to satisfy the #minimum will never cluster effectively
cluster_data <- cluster_data[apply(cluster_data, 2, sum, drop = F) >= min_core_genes,]

# Calculate distance matrix
dist_matrix <- dist(as.data.frame(t(cluster_data)))

# Hierarchical clustering
hc <- hclust(dist_matrix)

# Make DF to store iteratively tested clusters
cluster_df <- data.frame(matrix(NA, nrow = ncol(cluster_data), ncol = ncol(cluster_data)+1))
colnames(cluster_df) <- c(colnames(cluster_data), "Satisfied_clusters__")
cluster_df$Satisfied_clusters__=0

for(i in 1:nrow(cluster_df)){
  # Iteratively cut the dendrogram to obtain #clusters = i
  num_clusters <- i
  clusters <- cutree(hc, k = num_clusters)
  
  # store cluster assignments as row in dataframe
  cluster_df[i, 1:nrow(cluster_df)] <- clusters
  
  # for each total number of clusters, check if individual groupings have enough shared genes
  for(j in 1:num_clusters) {
    sample_subset <- cluster_data[, clusters == j, drop = F]
    length_subset <- ncol(sample_subset)
    n_samples_per_gene <- apply(sample_subset, 1, sum)
    
    num_satisfying <- sum(n_samples_per_gene == length_subset)
    if(num_satisfying >= min_core_genes) {
      cluster_df[i, 'Satisfied_clusters__'] <- cluster_df[i, 'Satisfied_clusters__'] + 1
    }
  }
}

# Only keep clusters with enough shared genes
cluster_df_filtered <- cluster_df[cluster_df$Satisfied_clusters__ == 1:nrow(cluster_df), ] 
cluster_df_filtered$Satisfied_clusters__ <- NULL

# Make list of sample names for each cluster (from the lowest acceptable cluster number)
sample_list <- lapply(1:max(cluster_df_filtered[1,]), function(i) {
    colnames(cluster_df_filtered)[cluster_df_filtered[1,] == i]
})

# Determine which cluster have both samples and references
has_ref_and_samp <- unlist(lapply(sample_list, function(ids) {
    any(ids %in% ref_ids) && any(ids %in% sample_ids) 
}))

# Make subset Pirate output for clusters with both samples and references
output_clusters <- lapply(sample_list[has_ref_and_samp], function(ids) {
    is_core_gene <- rowSums(cluster_data[, ids]) == length(ids)
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
    removed_ids <- unlist(sample_list[! has_ref_and_samp])
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
for (cluster in output_clusters) {
    out_dir_path <- file.path(args$fasta_output_path, paste0('cluster_', index))
    dir.create(out_dir_path, showWarnings = FALSE)
    for (gene_id in cluster$gene_family) {
        in_path <- file.path(args$gene_seq_dir_path, paste0(gene_id, '.fasta'))
        out_path <- file.path(out_dir_path, paste0(gene_id, '.fasta'))
        seqs <- read_fasta(in_path)
        seqs <- seqs[c(sample_ids, ref_ids)]
        seqs <- seqs[!is.na(seqs)]
        write_fasta(seqs, out_path)
        paste(seqs, out_path)
    }
}


