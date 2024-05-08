#!/usr/bin/env Rscript


# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/68/ea901183107b8c018617efda5b9fa7/Brady_comp.csv',
#     '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/tmp/09/5d1306ac2147a3b26630de2037d38b/Brady.csv',
#     'deleteme.csv',
#     '0.9'
# )
args <- as.list(args)
names(args) <- c("ani_matrix", "samp_ref_pairs", "out_path", "start_min_ani")
ani_matrix <- read.csv(args$ani_matrix, check.names = FALSE)
rownames(ani_matrix) <- as.character(colnames(ani_matrix)) #
samp_ref_pairs <- read.csv(args$samp_ref_pairs, header = FALSE, col.names = c("sample_id", "reference"))
samp_ref_pairs$sample_id <- as.character(samp_ref_pairs$sample_id) #
start_min_ani <- as.numeric(args$start_min_ani) # The minimum ANI for a reference to be assigned to a samples
end_min_ani <- max(c(0, start_min_ani - 0.3)) # How low the minimum can go if no samples can be assigned
ani_interval <- 0.05 # How much the minimum ANI threshold changes each time it is decreased

if (all(! is.na(samp_ref_pairs$reference))) {
    output <- samp_ref_pairs
} else {
    # Classify samples and references
    sample_ids <- samp_ref_pairs$sample_id
    user_ref_ids <- unique(samp_ref_pairs$reference[! is.na(samp_ref_pairs$reference)])
    down_ref_ids <- colnames(ani_matrix)[! colnames(ani_matrix) %in% c(sample_ids, user_ref_ids)]

    # Subset ANI matrix for samples without user-defined references
    refless_sample_ids <- samp_ref_pairs$sample_id[is.na(samp_ref_pairs$reference)]
    sample_ani <- ani_matrix[refless_sample_ids, c(user_ref_ids, down_ref_ids)]

    # Initialize empty data frame to hold the final assignments
    output <- data.frame(sample_id = character(0), reference = character(0))

    # Function to assign references given a minimum ANI
    assign_ref <- function(min_ani, table) {
        output <- NULL  # Reinitialize output here
        tf_table = table >= min_ani
        get_ref_max_samp <- function(ids) {
            n_samps <- apply(tf_table[, ids, drop = FALSE], MARGIN = 2, sum, na.rm = TRUE)
            if (max(n_samps) > 0) {
                best_ref <- names(n_samps)[which.max(n_samps)]
                samples_assigned <- rownames(tf_table[tf_table[, best_ref], best_ref, drop = FALSE])
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
    min_ani <- start_min_ani
    while (nrow(sample_ani) > 0 & min_ani >= end_min_ani) {
        # print(paste("Number of remaining samples: ", nrow(sample_ani)))
        rows_to_add <- assign_ref(min_ani, sample_ani)
        if (!is.null(rows_to_add) && nrow(rows_to_add) > 0) {
            sample_ani <- sample_ani[! rownames(sample_ani) %in% rows_to_add$sample_id, ]
            output <- rbind(output, rows_to_add)
        } else {
            min_ani <- min_ani - ani_interval
            # print(paste("Reducing min_ani to: ", min_ani))
        }
    }

    # If no reference found, fill with NAs
    if (nrow(output) == 0) {
        output <- data.frame(sample_id = refless_sample_ids, reference = rep(NA_character_, length(refless_sample_ids)))
    }

    # Combine with samples having user-defined references
    output <- rbind(output, samp_ref_pairs[!is.na(samp_ref_pairs$reference), ])
    rownames(output) <- output$sample_id

    # Sort output to match order of input samples
    output <- output[samp_ref_pairs$sample_id, ]

    # Convert NAs to a string not likely to be confused with a real ref ID
    output$reference[is.na(output$reference)] <- "__NULL__"
}

# Write output to file
write.table(output, file = args$out_path, sep = ',', quote = FALSE, row.names = FALSE, col.names = FALSE)
