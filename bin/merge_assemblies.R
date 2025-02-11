#!/usr/bin/env -S Rscript --vanilla

# Parse input TSVs
args <- commandArgs(trailingOnly = TRUE)
args <- c(
    '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/e4/bfab12294c9c2d87b08658aaa53302/Peronosporaceae.tsv',
    '/media/fosterz/external_primary/files/projects/work/current/pathogensurveillance/work/be/45c978e3382dfea3756890ced44802/Pythiaceae.tsv'
)
tsv_paths <- as.list(args)
tsv_data <- lapply(tsv_paths, read.csv, sep = '\t')

# Combine data for each taxon into a single table
out_data <- do.call(rbind, tsv_data)

# Add column for modified ID
modified_id <- gsub(out_data$LastMajorReleaseAccession, pattern = '[\\/:*?"<>| .]', replacement = '_')
out_data <- cbind(reference_id = modified_id, out_data)

# Write output table
write.table(out_data, file = 'merged_assembly_stats.tsv', sep = '\t', row.names = FALSE)

