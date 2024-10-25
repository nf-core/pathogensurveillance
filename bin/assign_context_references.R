#!/usr/bin/env Rscript

# This script attempts to find a minimal subset of references that have a range of similarity to each sample.
# This is done by making a list (data.frame) of "bins" representing a range of scaled ANI values for each sample,
# calculating which references can satisfy each bin, and selecting the references that satisfy the most bins
# until all bins are satisfied.

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     '/home/fosterz/projects/pathogensurveillance/work/db/03da7a3dc6394378a8cdc86e28060f/subgroup_comp.csv',
#     '/home/fosterz/projects/pathogensurveillance/work/db/03da7a3dc6394378a8cdc86e28060f/subgroup.csv',
#     '3',
#     '5',
#     'all_context_refs.csv'
# )
names(args) <- c('ani_matrix', 'sample_data', 'n_refs_closest', 'n_refs_contextual', 'output_path')
args <- as.list(args)
ani_matrix <- read.csv(args$ani_matrix, header = TRUE, check.names = FALSE)
rownames(ani_matrix) <- colnames(ani_matrix)
n_refs_closest <- as.integer(args$n_refs_closest)
n_refs_contextual <- as.integer(args$n_refs_contextual)

# Read sample data with user-defined references
sample_data <- read.csv(args$sample_data, header = FALSE, col.names = c('sample_id', 'references', 'usage'))
sample_ids <- unique(sample_data$sample_id)

# If 'exclusive' references are present for a sample, remove all other references. Also remove 'excluded' references
sample_data <- do.call(rbind, lapply(split(sample_data, sample_data$sample_id), function(table) {
    if (any(table$usage == 'exclusive')) {
        table <- table[table$usage == 'exclusive', , drop = FALSE]
    }
    excluded_ids <- table$sample_id[table$usage == 'excluded']
    table <- table[! table$usage %in% excluded_ids, , drop = FALSE]
    return(table)
}))
rownames(sample_data) <- NULL
all_ids <- unique(c(sample_data$references, sample_data$sample_id))
ani_matrix <- ani_matrix[row.names(ani_matrix) %in% all_ids, names(ani_matrix) %in% all_ids, drop = FALSE]

# Scale ANI values for each sample
rescale <- function(x) {
    if (length(x) > 1) {
        return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
    } else {
        return(x)
    }
}
ref_ids <- unique(rownames(ani_matrix)[! rownames(ani_matrix) %in% sample_ids])
ani_ref_v_samples <- ani_matrix[ref_ids, sample_ids, drop = FALSE]
ani_ref_v_samples[ani_ref_v_samples == 0] <- NA
ani_scaled <- as.data.frame(apply(ani_ref_v_samples, MARGIN = 2, rescale, simplify = FALSE), check.names = FALSE)

# Initialize list of selected references with closest references and required references
closest_refs <- unlist(lapply(sample_ids, function(id) {
    rownames(ani_scaled)[tail(order(ani_scaled[, id]), n = n_refs_closest)]
}))
required_refs <- unlist(lapply(sample_ids, function(id) {
    sample_data$references[sample_data$sample_id == id & sample_data$usage %in% c('required', 'exclusive')]
}))
selected_refs <- unique(c(closest_refs, required_refs))

# Initialized data.frame of bins that need a reference
bins <- seq(from = 0, to = 1, length.out = n_refs_contextual + 1)
bin_data <- data.frame(
    sample_id = rep(sample_ids, each = n_refs_contextual),
    from = rep(bins[1:(length(bins) - 1)], length(sample_ids)),
    to = rep(bins[2:length(bins)], length(sample_ids))
)
bin_data$refs <- lapply(1:nrow(bin_data), function(i) {
    is_in_bin <- ani_scaled[, bin_data$sample_id[i]] > bin_data$from[i] & ani_scaled[, bin_data$sample_id[i]] <= bin_data$to[i]
    rownames(ani_scaled)[! is.na(is_in_bin) & is_in_bin]
})

# Remove bins that have no references that can work
bin_data <- bin_data[unlist(lapply(bin_data$refs, length)) != 0, ]

# Remove bins that have references already chosen as "closest" references
filter_bins <- function(data, selected) {
    has_ref <- unlist(lapply(data$refs, function(refs) {
        any(refs %in% selected)
    }))
    data[! has_ref, ]
}
bin_data <- filter_bins(bin_data, selected_refs)

# Select reference that works for the most remaining bins and repeat until no bins are left
while (nrow(bin_data) > 0) {
    bin_counts <- table(unlist(bin_data$refs))
    best_ref <- names(bin_counts)[which.max(bin_counts)]
    selected_refs <- c(selected_refs, best_ref)
    bin_data <- filter_bins(bin_data, selected_refs)
}

# Write selected references to output file with one reference ID per line
writeLines(selected_refs, args$output_path)
