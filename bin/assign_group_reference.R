#!/usr/bin/env -S Rscript --vanilla

# Options
min_ani <- 90

# Parse inputs
# args <- commandArgs(trailingOnly = TRUE)

args <- c("/media/fosterz/external_primary/files/projects/work/current/nf-core-plantpathsurveil/work/ab/6cedf39ed4ed7155dce19a44ac6b58/comp.csv",
          '/media/fosterz/external_primary/files/projects/work/current/nf-core-plantpathsurveil/work/tmp/05/aa27c0075d6f771079cf3f8af9e954/"xan_test.csv',
          "result.tsv")
args <- as.list(args)
names(args) <- c("ani_matrix", "samp_ref_pairs", "out_path")
ani_matrix <- read.csv(args$ani_matrix)
samp_ref_pairs <- read.csv(args$samp_ref_pairs)

# Convert triangular distance matrix to sample vs reference matrix

# Remove samples with references already assigned by the user

# Create function to assign as many samples to a reference as possible

# Repeat execution of function until no samples are left

# Put samples back into order they came in


# Save which samples have which references in an output CSV
write.table(samp_ref_pairs, file = args$out_path, sep = '\t', quote = FALSE, row.names = FALSE)
