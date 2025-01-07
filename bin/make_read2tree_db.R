#!/usr/bin/env Rscript

# This script converts the output of BUSCO on reference sequences into a reference database for use with read2tree.
# Since read2tree requires references to be labeled with 5 letter species codes, this script will produce arbitrary
# 5 letter codes as well as a lookup table to correlate codes with their reference ID.

# Input format:
#    * name of lineage DB used for BUSCO (e.g. 'eukaryota_odb10')
#    * name of DNA DB output file
#    * name of AA DB output folder
#    * name of ref metadata output file
#    * one or more busco output directories

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     'eukaryota_odb10',
#     'prefix_dna_ref.fa',
#     'prefix_r2t_markers',
#     'prefix_rt2_ref_data.csv',
#     '/home/fosterz/data/files/projects/current/pathogensurveillance/work/5d/d22a2f4b743c17b97b7f0e4f920133/GCA_019155715_1-eukaryota_odb10-busco',
#     '/home/fosterz/data/files/projects/current/pathogensurveillance/work/5d/d22a2f4b743c17b97b7f0e4f920133/GCA_000365505_1-eukaryota_odb10-busco',
#     '/home/fosterz/data/files/projects/current/pathogensurveillance/work/5d/d22a2f4b743c17b97b7f0e4f920133/GCA_020800235_1-eukaryota_odb10-busco'
# )
args <- as.list(args)
lineage_db <- args[[1]]
dna_ref_out_path <- args[[2]]
aa_ref_out_path <- args[[3]]
ref_data_out_path <- args[[4]]
busco_dir_paths <- unlist(args[5:length(args)])

# Assign random 5-letter codes for references
number_to_letter_code <- function(number, base = 26, suffix = "") {
    # From https://stackoverflow.com/questions/44269918/numeric-to-alphabetic-lettering-function-in-r
    number1 <- number - 1
    last_digit <- number1 %% base
    rest <- number1 %/% base
    suffix <- paste0(LETTERS[last_digit + 1], suffix)
    if (rest > 0) Recall(rest, base, suffix) else suffix
}
pad_strings <- function(x, len = 5, pad = "A") {
    to_add <- unlist(lapply(len - nchar(x), function(n_pad) paste0(rep(pad, n_pad), collapse = '')))
    paste0(to_add, x)
}
ref_data <- data.frame(
    ref_id = sub(basename(busco_dir_paths), pattern = paste0('-', lineage_db, '-busco'), replacement = ''),
    r2t_ref_id = pad_strings(unlist(lapply(seq_along(busco_dir_paths), number_to_letter_code))),
    busco_dir_path = busco_dir_paths
)

# Read and reformat FASTA data
read_and_reformat_fasta <- function(path_pattern) {
    do.call(rbind, lapply(1:nrow(ref_data), function(i) {
        fasta_paths <- Sys.glob(file.path(ref_data$busco_dir_path[i], path_pattern))
        do.call(rbind, lapply(fasta_paths, function(fasta_path) {
            lines <- readLines(fasta_path)
            busco_id <- toupper(tools::file_path_sans_ext(basename(fasta_path)))
            gene_id <- paste0(ref_data$r2t_ref_id[i], '-', busco_id)
            ref_id <- ref_data$ref_id[i]
            lines[1] <- paste0('>', gene_id, ' | ', busco_id, ' | ', gene_id, ' | ', paste0('[', ref_id, ']'))
            data.frame(
                busco_id = busco_id,
                gene_id = gene_id,
                ref_id = ref_id,
                fasta = paste0(lines, collapse = '\n')
            )
        }))
    }))
}
single_copy_fna_data <- read_and_reformat_fasta("*/run_*/busco_sequences/single_copy_busco_sequences/*.fna")
single_copy_faa_data <- read_and_reformat_fasta("*/run_*/busco_sequences/single_copy_busco_sequences/*.faa")

# Save FASTA file with all DNA sequences
writeLines(single_copy_fna_data$fasta, dna_ref_out_path)

# Save FASTA files of AA sequences grouped into files by locus
dir.create(aa_ref_out_path, showWarnings = FALSE)
for (gene_data in split(single_copy_faa_data, single_copy_faa_data$busco_id)) {
    writeLines(gene_data$fasta, file.path(aa_ref_out_path, paste0(gene_data$busco_id[1], '.fasta')))
}

# Save table with info on which random code corresponds to which reference ID
write.table(ref_data, file = ref_data_out_path, sep = ',', quote = FALSE, row.names = FALSE)

