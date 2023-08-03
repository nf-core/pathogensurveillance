#!/usr/bin/env -S Rscript --vanilla

# Options
min_ani <- 0.9

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c("/media/fosterz/external_primary/files/projects/work/current/nf-core-plantpathsurveil/work/e1/4d319311c6e90eacb50cbffc7b9932/comp.csv",
#           '/media/fosterz/external_primary/files/projects/work/current/nf-core-plantpathsurveil/work/tmp/81/91384f8b280e1ae1078a76c99934ab/subgroup.csv',
#           "result.tsv")
args <- as.list(args)
names(args) <- c("ani_matrix", "samp_ref_pairs", "out_path")
ani_matrix <- read.csv(args$ani_matrix, check.names = FALSE)
rownames(ani_matrix) <- colnames(ani_matrix)
samp_ref_pairs <- read.csv(args$samp_ref_pairs, header = FALSE, col.names = c("sample_id", "reference"))

# Infer IDs for samples, user-defined references, and downloaded references
sample_ids <- samp_ref_pairs$sample_id
user_ref_ids <- unique(samp_ref_pairs$reference[! is.na(samp_ref_pairs$reference)])
down_ref_ids <- colnames(ani_matrix)[! colnames(ani_matrix) %in% c(sample_ids, user_ref_ids)]
ref_ids <- c(user_ref_ids, down_ref_ids)

# Convert triangular distance matrix to sample vs reference matrix
refless_sample_ids <- samp_ref_pairs$sample_id[is.na(samp_ref_pairs$reference)]
sample_ani <- ani_matrix[refless_sample_ids, ref_ids]

# Create function to assign as many samples to a reference as possible
assign_ref <- function(table) {
  tf_table = table >= min_ani
  
  get_ref_max_samp <- function(ids) {
    n_samps <- apply(tf_table[, ids, drop = FALSE], MARGIN = 2, sum, simplify = TRUE)
    if (max(n_samps) > 0) {
      best_ref <- names(n_samps)[which.max(n_samps)]
      samples_assigned <- names(tf_table[tf_table[, best_ref], best_ref])
      return(data.frame(sample_id = samples_assigned, reference = best_ref))
    }
  }
  
  # Check for best user-defined reference
  output = NULL
  if (length(user_ref_ids) > 0) {
    output <- get_ref_max_samp(user_ref_ids)
  }
  
  # Check for best downloaded reference
  if (is.null(output)) {
    output <- get_ref_max_samp(down_ref_ids)
  }
  
  return(output)
}

# Repeat execution of function until no samples are left
output <- list()
while (nrow(sample_ani) > 0) {
  rows_to_add <- assign_ref(sample_ani)
  sample_ani <- sample_ani[! rownames(sample_ani) %in% rows_to_add$sample_id, ]
  output <- c(output, list(rows_to_add))
}
output <- do.call(rbind, output)

# Add back samples with user-defined references
output <- rbind(output, samp_ref_pairs[!is.na(samp_ref_pairs$reference), ])
rownames(output) <- output$sample_id

# Put samples back into order they came in
output <- output[samp_ref_pairs$sample_id, ]

# Save which samples have which references in an output CSV
write.table(output, file = args$out_path, sep = ',', quote = FALSE, row.names = FALSE)
