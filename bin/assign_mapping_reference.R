#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     '/home/fosterz/projects/pathogensurveillance/work/6f/a8638157f5028f753f9ff512bb3de0/subgroup_comp.csv',
#     '/home/fosterz/projects/pathogensurveillance/work/6f/a8638157f5028f753f9ff512bb3de0/subgroup.csv',
#     'deleteme.csv',
#     '0.95'
# )
args <- as.list(args)
names(args) <- c("ani_matrix", "samp_ref_pairs", "out_path", "start_min_ani")
ani_matrix <- read.csv(args$ani_matrix, check.names = FALSE)
rownames(ani_matrix) <- as.character(colnames(ani_matrix))
samp_ref_pairs <- read.csv(args$samp_ref_pairs, header = FALSE, col.names = c("sample_id", "ref_id", "ref_name", "ref_desc", "usage"), sep = '\t')
samp_ref_pairs$sample_id <- as.character(samp_ref_pairs$sample_id)
start_min_ani <- as.numeric(args$start_min_ani) # The minimum ANI for a reference to be assigned to a samples
end_min_ani <- max(c(0, start_min_ani - 0.3)) # How low the minimum can go if no samples can be assigned
ani_interval <- 0.05 # How much the minimum ANI threshold changes each time it is decreased

# If 'exclusive'/ 'required' references are present for a sample, remove all other references. Also remove 'excluded' references
samp_ref_pairs <- do.call(rbind, lapply(split(samp_ref_pairs, samp_ref_pairs$sample_id), function(sample_data) {
    if (any(sample_data$usage == 'exclusive')) {
        sample_data <- sample_data[sample_data$usage == 'exclusive', , drop = FALSE]
    }
    if (any(sample_data$usage == 'required')) {
        sample_data <- sample_data[sample_data$usage == 'required', , drop = FALSE]
    }
    excluded_ids <- sample_data$sample_id[sample_data$usage == 'excluded']
    sample_data <- sample_data[! sample_data$usage %in% excluded_ids, , drop = FALSE]
    return(sample_data)
}))
rownames(samp_ref_pairs) <- NULL

# Function to assign references given a minimum ANI
reference_ids <- unique(samp_ref_pairs$ref_id)
assign_ref <- function(sample_ids, min_ani) {
    valid_samples_for_ref <- lapply(reference_ids, function(ref_id) {
        good_ani <- ani_matrix[ref_id, sample_ids] >= min_ani
        can_use_ref <- unlist(lapply(sample_ids, function(sample_id) {
            any(samp_ref_pairs$sample_id == sample_id & samp_ref_pairs$ref_id == ref_id)
        }))
        sample_ids[good_ani & can_use_ref]
    })
    names(valid_samples_for_ref) <- reference_ids
    n_samples <- unlist(lapply(valid_samples_for_ref, length))
    if (max(n_samples) == 0) {
        return(NULL)
    }
    best_refs <- reference_ids[n_samples == max(n_samples)]
    mean_ani <- vapply(best_refs, function(ref_id) {
        mean(unlist(ani_matrix[ref_id, valid_samples_for_ref[[ref_id]], drop = FALSE]))
    }, FUN.VALUE = numeric(1))
    best_ref <- best_refs[which.max(mean_ani)]
    return(data.frame(sample_id = valid_samples_for_ref[[best_ref]], reference_id = best_ref))
}

# Main loop for assignment based on given min_ani
unassigned_sample_ids <- unique(samp_ref_pairs$sample_id)
output <- data.frame(sample_id = character(0), reference_id = character(0))
min_ani <- start_min_ani
while (length(unassigned_sample_ids) > 0 & min_ani >= end_min_ani) {
    # print(paste("Number of remaining samples: ", length(unassigned_sample_ids)))
    rows_to_add <- assign_ref(unassigned_sample_ids, min_ani)
    if (! is.null(rows_to_add) && nrow(rows_to_add) > 0) {
        unassigned_sample_ids <- unassigned_sample_ids[! unassigned_sample_ids %in% rows_to_add$sample_id]
        output <- rbind(output, rows_to_add)
    } else {
        min_ani <- min_ani - ani_interval
        # print(paste("Reducing min_ani to: ", min_ani))
    }
}

# If no reference found, fill with NAs
if (length(unassigned_sample_ids) > 0) {
    output <- rbind(output, data.frame(sample_id = unassigned_sample_ids, reference_id = "__NULL__"))
}

# Write output to file
write.table(output, file = args$out_path, sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)
