#!/usr/bin/env Rscript

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

# Convert BUSCO output directories into the unified gene_sequences/ format.
# Args: <output_dir> <busco_dir_1> [<busco_dir_2> ...]
# Output: <output_dir>/gene_sequences/<busco_id>.fasta
#   Headers are genome IDs; duplicate headers indicate multi-copy genes.
# Fragmented BUSCO sequences are excluded.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: format_busco_genes.R <output_dir> <busco_dir_1> [<busco_dir_2> ...]")
}
out_dir <- args[1]
busco_dirs <- args[2:length(args)]

# Infer sample IDs from BUSCO directory names
run_paths <- list.files(busco_dirs[1], recursive = TRUE, include.dirs = TRUE, pattern = '^run_')
if (length(run_paths) == 0) {
    stop("No run_ directory found in ", busco_dirs[1])
}
lineage_db_id <- sub("^run_", "", basename(run_paths[1]))
names(busco_dirs) <- sub(paste0('--', lineage_db_id, '-busco$'), "", basename(busco_dirs))

seq_dir <- file.path(out_dir, "gene_sequences")
dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

read_fasta <- function(file_path) {
    lines <- readLines(file_path)
    lines <- lines[lines != ""]
    seqs <- list()
    current_name <- NULL
    current_buffer <- character(0)
    for (line in lines) {
        if (grepl("^>", line)) {
            if (!is.null(current_name)) {
                seqs[[current_name]] <- paste(current_buffer, collapse = "")
            }
            current_name <- sub("^>([^ ]+).*", "\\1", line)
            current_buffer <- character(0)
        } else {
            current_buffer <- c(current_buffer, toupper(line))
        }
    }
    if (!is.null(current_name)) {
        seqs[[current_name]] <- paste(current_buffer, collapse = "")
    }
    return(seqs)
}

for (i in seq_along(busco_dirs)) {
    sample_id <- names(busco_dirs)[i]
    busco_dir <- busco_dirs[i]

    # Find all BUSCO gene FASTA files recursively
    fasta_files <- list.files(busco_dir, recursive = TRUE, pattern = "\\.fna$", full.names = TRUE)

    for (fasta_file in fasta_files) {
        # Include single_copy and multi_copy; exclude fragmented and other files
        if (grepl("single_copy_busco_sequences", fasta_file)) {
            # keep
        } else if (grepl("multi_copy_busco_sequences", fasta_file)) {
            # keep
        } else {
            next
        }

        gene_id <- sub("\\.fna$", "", basename(fasta_file))
        seqs <- read_fasta(fasta_file)
        out_path <- file.path(seq_dir, paste0(gene_id, ".fasta"))
        h <- file(out_path, open = "at")
        for (seq_name in names(seqs)) {
            writeLines(paste0(">", sample_id), h)
            writeLines(seqs[[seq_name]], h)
        }
        close(h)
    }
}

cat("format_busco_genes:\n  r-base: ", R.version.string, "\n", sep = "")
