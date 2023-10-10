#!/usr/bin/env -S Rscript --vanilla

start_min_ani <- 0.9 # The minimum ANI for a reference to be assigned to a samples 
end_min_ani <- 0.7 # How low the minimum can go if no samples can be assigned
ani_interval <- 0.02 # How much the minimum ANI threshold changes each time it is decreased 

# Function to assign references based on a given min_ani threshold
assign_ref_based_on_min_ani <- function(min_ani, sample_ani, user_ref_ids, down_ref_ids) {
  # Initialize empty data frame to hold the final assignments
  output <- data.frame(sample_id = character(0), reference = character(0))
  
  # Internal function to assign references
  assign_ref <- function(min_ani, table) {
    output <- NULL  # Reinitialize output here
    tf_table = table >= min_ani
    
    get_ref_max_samp <- function(ids) {
      n_samps <- apply(tf_table[, ids, drop = FALSE], MARGIN = 2, sum, na.rm = TRUE)
      if (max(n_samps) > 0) {
        best_ref <- names(n_samps)[which.max(n_samps)]
        samples_assigned <- names(tf_table[tf_table[, best_ref], best_ref])
        if (length(samples_assigned) > 0) {
          return(data.frame(sample_id = samples_assigned, reference = best_ref))
        }
      }
      return(NULL)  # Return NULL if no valid data.frame can be created
    }
    
    if (length(user_ref_ids) > 0) {
      output <- get_ref_max_samp(user_ref_ids)
    }
    
    if (is.null(output)) {
      output <- get_ref_max_samp(down_ref_ids)
    }
    
    return(output)
  }
  
  # Main loop for assignment based on given min_ani
  while (nrow(sample_ani) > 0) {
    print(paste("Number of remaining samples: ", nrow(sample_ani)))  # New print statement
    rows_to_add <- assign_ref(min_ani, sample_ani)
    
    if (!is.null(rows_to_add) && nrow(rows_to_add) > 0) {
      sample_ani <- sample_ani[! rownames(sample_ani) %in% rows_to_add$sample_id, ]
      output <- rbind(output, rows_to_add)
    } else {
      break
    }
  }
  
  return(output)
}

# Parse inputs

args <- commandArgs(trailingOnly = TRUE)
args <- as.list(args)
names(args) <- c("ani_matrix", "samp_ref_pairs", "out_path")
ani_matrix <- read.csv(args$ani_matrix, check.names = FALSE)
rownames(ani_matrix) <- colnames(ani_matrix)
samp_ref_pairs <- read.csv(args$samp_ref_pairs, header = FALSE, col.names = c("sample_id", "reference"))

# Classify samples and references
sample_ids <- samp_ref_pairs$sample_id
user_ref_ids <- unique(samp_ref_pairs$reference[! is.na(samp_ref_pairs$reference)])
down_ref_ids <- colnames(ani_matrix)[! colnames(ani_matrix) %in% c(sample_ids, user_ref_ids)]

# Subset ANI matrix for samples without user-defined references
refless_sample_ids <- samp_ref_pairs$sample_id[is.na(samp_ref_pairs$reference)]
sample_ani <- ani_matrix[refless_sample_ids, c(user_ref_ids, down_ref_ids)]

# Initialize empty data frame for final output
final_output <- data.frame(sample_id = character(0), reference = character(0))

# Main loop: loop over min_ani values from 0.9 to 0.7
if (length(refless_sample_ids) > 0) {
    min_ani <- start_min_ani
    while (min_ani >= end_min_ani) {
        print(paste("Current min_ani: ", min_ani))  # Debugging print statement
        temp_output <- assign_ref_based_on_min_ani(min_ani, sample_ani, user_ref_ids, down_ref_ids)
        if (nrow(temp_output) > 0) {
            final_output <- rbind(final_output, temp_output)
            break
        }
        min_ani <- min_ani - ani_interval
    }
}

# If no reference found, fill with NAs
if (nrow(final_output) == 0) {
  final_output <- data.frame(sample_id = refless_sample_ids, reference = rep(NA_character_, length(refless_sample_ids)))
}

# Combine with samples having user-defined references
final_output <- rbind(final_output, samp_ref_pairs[!is.na(samp_ref_pairs$reference), ])
rownames(final_output) <- final_output$sample_id

# Sort output to match order of input samples
final_output <- final_output[samp_ref_pairs$sample_id, ]

# Convert NAs to a string not likely to be confused with a real ref ID
final_output$reference[is.na(final_output$reference)] <- "__NULL__"

# Write output to file
write.table(final_output, file = args$out_path, sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)