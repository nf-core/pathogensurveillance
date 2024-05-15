#!/usr/bin/env Rscript

# This script attempts to find a minimal subset of references that have a range of similarity to each sample.
# This is done by making a list (data.frame) of "bins" representing a range of scaled ANI values for each sample,
# calculating which references can satisfy each bin, and selecting the references that satisfy the most bins
# until all bins are satisfied.

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- c('test/output/mycobacteroides_small/sourmash_compare/all_comp.csv', '...', '3', '5')
names(args) <- c('ani_matrix', 'metadata', 'n_refs_closest', 'n_refs_contextual')
args <- as.list(args)
ani_matrix <- read.csv(args$ani_matrix, header = TRUE, check.names = FALSE)
rownames(ani_matrix) <- colnames(ani_matrix)
n_refs_closest <- as.integer(args$n_refs_closest)
n_refs_contextual <- as.integer(args$n_refs_contextual)


sample_ids <- rownames(ani_matrix)[1:3] # TODO: base on input files
ref_ids <- rownames(ani_matrix)[4:23]

# Scale ANI values for each sample
rescale <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}
ani_scaled <- apply(ani_matrix[ref_ids, sample_ids], 2, rescale)

# Initialize list of selected references with closest references
selected_refs <- unique(unlist(lapply(sample_ids, function(id) {
    rownames(ani_scaled)[tail(order(ani_scaled[, id]), n = n_refs_closest)]
})))

# Initialized data.frame of bins that need a reference
bins <- seq(from = 0, to = 1, length.out = n_refs_contextual + 1)
bin_data <- data.frame(
    sample_id = rep(sample_ids, each = n_refs_contextual),
    from = rep(bins[1:(length(bins) - 1)], length(sample_ids)),
    to = rep(bins[2:length(bins)], length(sample_ids))
)
bin_data$refs <- lapply(1:nrow(bin_data), function(i) {
    is_in_bin <- between(ani_scaled[, bin_data$sample_id[i]], bin_data$from[i], bin_data$to[i])
    rownames(ani_scaled)[is_in_bin]
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

