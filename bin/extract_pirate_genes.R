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

# Extract per-gene-family FASTA sequences from PIRATE results and GFF files.
# Args: <pirate_results_dir> <output_dir>
# Output: <output_dir>/gene_sequences/<gene_family>.fasta
#   Headers are genome IDs; duplicate headers indicate multi-copy genes.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
    stop("Usage: extract_pirate_genes.R <pirate_results_dir> <output_dir>")
}
pirate_dir <- args[1]
out_dir <- args[2]

# pirate_dir <- "/home/fosterz/projects/pathogensurveillance/pathogensurveillance/work/b8/3b8c80853dab7019e16b84d027f00f/pirate_results"
# out_dir <- "deleteme"

# ---------------------------------------------------------------------------
# Parse PIRATE gene families table
# ---------------------------------------------------------------------------
pirate_file <- file.path(pirate_dir, "PIRATE.gene_families.ordered.tsv")
if (!file.exists(pirate_file)) {
    stop("PIRATE table not found: ", pirate_file)
}

pirate_lines <- readLines(pirate_file)
header_line <- pirate_lines[grepl("^allele_name", pirate_lines)]
if (length(header_line) == 0) {
    stop("No header line found in PIRATE table.")
}
header_parts <- strsplit(header_line, "\t")[[1]]

# Determine first genome column index (0-based)
idx <- 19
if (any(grepl("\tno_loci\t", header_line))) {
    idx <- 20
}
if (any(grepl("\tcluster_order\t", header_line))) {
    idx <- 22
}
genome_cols <- seq(idx + 1, length(header_parts)) # 1-based R indices
genome_ids <- header_parts[genome_cols]

# Build locus -> family and locus -> genome lookups
locus_to_family <- character(0)
locus_to_genome <- character(0)
for (line in pirate_lines[-1]) {
    if (line == "" || grepl("^#", line)) next
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) < max(genome_cols)) next
    family <- parts[2]
    for (g_i in seq_along(genome_cols)) {
        genome_id <- genome_ids[g_i]
        cell <- parts[genome_cols[g_i]]
        if (cell == "") next
        # Split on ; : /
        loci <- unlist(strsplit(cell, "[;:/]"))
        loci <- loci[loci != ""]
        for (locus in loci) {
            locus_to_family[locus] <- family
            locus_to_genome[locus] <- genome_id
        }
    }
}

families <- unique(locus_to_family)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
reverse_complement <- function(seq) {
    seq <- chartr("ACGTacgt", "TGCAtgca", seq)
    seq <- paste(rev(strsplit(seq, "")[[1]]), collapse = "")
    return(seq)
}

# ---------------------------------------------------------------------------
# Process each genome GFF and write per-family FASTA files
# ---------------------------------------------------------------------------
gff_dir <- file.path(pirate_dir, "modified_gffs")
if (!dir.exists(gff_dir)) {
    stop("modified_gffs directory not found: ", gff_dir)
}

gff_files <- list.files(gff_dir, pattern = "\\.gff$", full.names = TRUE)
if (length(gff_files) == 0) {
    stop("No GFF files found in ", gff_dir)
}

seq_dir <- file.path(out_dir, "gene_sequences")
dir.create(seq_dir, recursive = TRUE, showWarnings = FALSE)

for (gff_path in gff_files) {
    genome_id <- sub("\\.gff$", "", basename(gff_path))

    lines <- readLines(gff_path)

    # Locate FASTA section
    fasta_idx <- which(lines == "##FASTA")
    if (length(fasta_idx) == 0) {
        warning("No ##FASTA section in ", gff_path)
        next
    }
    fasta_idx <- fasta_idx[1]

    # Parse sequences
    seq_lines <- lines[(fasta_idx + 1):length(lines)]
    seq_lines <- seq_lines[seq_lines != ""]
    contig_seqs <- list()
    current_contig <- NULL
    current_buffer <- character(0)
    for (sline in seq_lines) {
        if (grepl("^>", sline)) {
            if (!is.null(current_contig)) {
                contig_seqs[[current_contig]] <- paste(current_buffer, collapse = "")
            }
            current_contig <- sub("^>([^ ]+).*", "\\1", sline)
            current_buffer <- character(0)
        } else {
            # Replace unexpected chars with N, uppercase
            sline <- toupper(sline)
            sline <- gsub("[^ACGTN-]", "N", sline)
            current_buffer <- c(current_buffer, sline)
        }
    }
    if (!is.null(current_contig)) {
        contig_seqs[[current_contig]] <- paste(current_buffer, collapse = "")
    }

    # Parse feature lines (before ##FASTA)
    feat_lines <- lines[1:(fasta_idx - 1)]
    for (fline in feat_lines) {
        if (fline == "" || grepl("^#", fline)) next
        parts <- strsplit(fline, "\t")[[1]]
        if (length(parts) < 9) next
        type <- parts[3]
        if (type == "gene") next

        contig <- parts[1]
        start <- as.integer(parts[4])
        end <- as.integer(parts[5])
        strand <- parts[7]

        attrs <- parts[9]
        attr_pairs <- strsplit(attrs, ";")[[1]]
        id <- NA_character_
        for (pair in attr_pairs) {
            if (grepl("^locus_tag=", pair)) {
                id <- sub("^locus_tag=", "", pair)
                break
            } else if (grepl("^ID=", pair)) {
                id <- sub("^ID=", "", pair)
                break
            }
        }
        if (is.na(id)) next
        if (!id %in% names(locus_to_family)) next

        family <- locus_to_family[id]

        if (!contig %in% names(contig_seqs)) {
            warning("Contig ", contig, " not found in FASTA section of ", gff_path)
            next
        }
        contig_seq <- contig_seqs[[contig]]
        seq <- substr(contig_seq, start, end)
        if (strand == "-") {
            seq <- reverse_complement(seq)
        }

        out_path <- file.path(seq_dir, paste0(family, ".fasta"))
        h <- file(out_path, open = "at")
        writeLines(paste0(">", genome_id), h)
        writeLines(seq, h)
        close(h)
    }
}

# Write versions info
cat("extract_pirate_genes:\n  r-base: ", R.version.string, "\n", sep = "")
