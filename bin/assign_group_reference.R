#!/usr/bin/env -S Rscript --vanilla

# Options
min_ani <- 90

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)

# args <- c("",
#           "",
#           "result.tsv")
args <- as.list(args)
names(args) <- c("ani_matrix", "samp_ref_pairs", "out_path")
ani_matrix <- read.csv(args$ani_matrix)
samp_ref_pairs <- read.csv(args$samp_ref_pairs)

# Save to output file
write.table(samp_ref_pairs, file = args$out_path, sep = '\t', quote = FALSE, row.names = FALSE)
