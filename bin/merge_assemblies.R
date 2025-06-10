#!/usr/bin/env -S Rscript --vanilla

# MIT License
#
# Copyright (c) Zachary S.L. Foster and Niklaus J. Grunwald
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


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
