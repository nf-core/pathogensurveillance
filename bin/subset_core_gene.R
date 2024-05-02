#!/usr/bin/env Rscript
setwd('/home/buchananr/Desktop/pathogensurveillance/subset_core_inputs_for_logan/')
# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- c('subgroup.tsv',
          'subgroup_feat_seqs_renamed',
          'samplesheet.valid.csv',
          '10', '0.8', '0.5', '100', 'subgroup_core_genes.tsv', 'subgroup_feat_seqs')
names(args) <- c("gene_families", "gene_seq_dir_path", "metadata", "min_core_genes", "min_core_samps", "min_core_refs", "max_core_genes", "csv_output_path", "fasta_output_path")
args <- as.list(args)
raw_gene_data <- read.csv(args$gene_families, header = TRUE, sep = '\t', check.names = FALSE)
metadata <- read.csv(args$metadata, header = TRUE, sep = ',', row.names = NULL, check.names = FALSE)
min_core_genes <- as.integer(args$min_core_genes)
max_core_genes <- as.integer(args$max_core_genes)
min_core_samps <- as.numeric(args$min_core_samps) #previously as.integer. was storing as 0
min_core_refs <- as.numeric(args$min_core_refs) #previously as.integer. was storing as 0

raw_gene_data1 <- raw_gene_data

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

gene_data_clust <- gene_data[,colnames(gene_data)%in%all_ids]

########
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
  print(paste0('samples:  ', length(current_sample_ids)))
  print(paste0('refs:     ', length(current_ref_ids)))
  print(paste0('core:     ', current_core_genes))
  
  # Find how many genes each sample are missing or multi-copy
  n_bad_samp <- colSums(gene_data_subset[, current_sample_ids, drop = FALSE] != 1)
  n_bad_ref <- colSums(gene_data_subset[, current_ref_ids, drop = FALSE] != 1)
  
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
#######
# Save the IDs of any samples and references that were removed from the analysis
removed_sample_ids <- sample_ids[! sample_ids %in% current_sample_ids]
removed_ref_ids <- ref_ids[! ref_ids %in% current_ref_ids]
writeLines(removed_sample_ids, con = 'removed_sample_ids.txt')
writeLines(removed_ref_ids, con = 'removed_ref_ids.txt')

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

passing_sample_ids
for (gene_id in output$gene_family) {
  in_path <- file.path(args$gene_seq_dir_path, paste0(gene_id, '.fasta'))
  out_path <- file.path(args$fasta_output_path, paste0(gene_id, '.fasta'))
  seqs <- read_fasta(in_path)
  seqs <- seqs[passing_sample_ids]
  seqs <- seqs[!is.na(seqs)]
  write_fasta(seqs, out_path)
  paste(seqs, out_path)
}









########
#Splitting trees if specified

cluster_data <- raw_gene_data1[,23:ncol(raw_gene_data1)]
cluster_data[, all_ids] <- lapply(cluster_data[, all_ids], function(column) {
  unlist(lapply(strsplit(column, split = '[;:]'), length))
})

# treating multiple copies as missing
cluster_data[cluster_data >1] <- 0

# Samples that don't have enough genes to satisfy the #minimum will never cluster effectively
cluster_data <- cluster_data[apply(cluster_data, 2, sum, drop=F) >= min_core_genes,]

# Calculate distance matrix
cluster_data1=as.data.frame(t(cluster_data))
dist_matrix <- dist(cluster_data1)

# Hierarchical clustering
hc <- hclust(dist_matrix)

# Make DF to store iteratively tested clusters
cluster_df <- data.frame(matrix(NA, nrow = ncol(cluster_data), ncol = ncol(cluster_data)+1))
colnames(cluster_df) <- c(colnames(cluster_data), "Satisfied_clusters")
cluster_df$Satisfied_clusters=0

for(i in 1:ncol(cluster_data)){
  # Iteratively cut the dendrogram to obtain #clusters = i
  num_clusters <- i
  clusters <- cutree(hc, k = num_clusters)
  
  # store cluster assignments as row in dataframe
  cluster_df[i,1:nrow(cluster_df)] <- clusters
  
  # for each total number of clusters, check if individual groupings have enough shared genes
  for(j in 1:num_clusters){
    row_subset <- cluster_data[which(clusters==j),,drop=F]
    length_subset <- nrow(row_subset)
    row_sums <- apply(row_subset, 2, sum)
    
    num_satisfying <- length(which(row_sums==length_subset))
    if(num_satisfying >= min_core_genes) {
      cluster_df[i,'Satisfied_clusters'] <- cluster_df[i,'Satisfied_clusters']+1
    }
  }
}

# only keep clusters with enough shared genes
cluster_df$Satisfied_clusters <- (cluster_df$Satisfied_clusters / 1:nrow(cluster_df)-1)
cluster_df1 <- cluster_df[cluster_df$Satisfied_clusters==0,] # 0 if satisfied
cluster_df1$Satisfied_clusters <- NULL

for(i in 1:nrow(cluster_df1)){
  
  for(j in 1:length(unique(as.character(cluster_df1[i,])))){
    
    # for #cluster = i, grab each cluster
    col_vars <- which(cluster_df1[i,] == j)
    
    # require at least one sample and one reference in order to keep cluster, otherwise convert values to NA
    if( (length (which (colnames (cluster_df1)[col_vars] %in% ref_ids))>1
         && length (which (colnames (cluster_df1)[col_vars] %in% sample_ids))>1)==F){
      cluster_df1[i,col_vars] <- NA
    }
  }
}

sample_list <- list()

# make list of sample names for each cluster (from the lowest acceptable cluster number)
for(i in 1:length(unique(as.character(cluster_df1[1,])))){
  sample_list[[i]] <- colnames(cluster_df1)[cluster_df1[1,]==i]
}

sample_list
