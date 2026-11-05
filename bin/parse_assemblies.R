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


library(RcppSimdJson)

# Parse taxonomy inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     '~/projects/pathogensurveillance/path_surveil_data/assembly_metadata/Carnobacteriaceae.json',
#     'deleteme.csv'
# )
args <- as.list(args)
in_path <- args[[1]]
out_path <- args[[2]]

# Add header to output file
header <- data.frame(
    reference_id = character(0),
    accession = character(0),
    assembly_level = character(0),
    assembly_status = character(0),
    assembly_type = character(0),
    hosts = character(0),
    organism_name = character(0),
    tax_id = character(0),
    contig_l50 = numeric(0),
    contig_n50 = numeric(0),
    coverage = numeric(0),
    number_of_component_sequences = numeric(0),
    number_of_contigs = numeric(0),
    total_ungapped_length = numeric(0),
    total_sequence_length = numeric(0),
    source_database = character(0),
    is_type = logical(0),
    is_annotated = logical(0),
    is_atypical = logical(0),
    checkm_completeness = numeric(0),
    checkm_contamination = numeric(0)
)
write.table(header, file = out_path, sep = '\t', quote = FALSE)

# Parse and write assembly metadata one at a time
in_handle <- file(in_path)
open(in_handle)
out_handle <- file(out_path, open = 'a')
while (length(line <- readLines(in_handle, n = 1)) != 0) {
    assem_data <- RcppSimdJson::fparse(line, always_list = TRUE)[[1]]
    attributes <- assem_data$assembly_info$biosample$attributes
    hosts <- paste0(attributes$value[attributes$name == 'host'], collapse = ';')
    data_parts <- list(
        reference_id = gsub(assem_data$accession, pattern = '[\\/:*?"<>| .]', replacement = '_'),
        accession = assem_data$accession,
        assembly_level = assem_data$assembly_info$assembly_level,
        assembly_status = assem_data$assembly_info$assembly_status,
        assembly_type = assem_data$assembly_info$assembly_type,
        hosts = ifelse(hosts == '', NA_character_, hosts),
        organism_name = gsub(assem_data$organism$organism_name, pattern = '\\[|\\]', replacement = ''),
        tax_id = as.character(assem_data$organism$tax_id),
        contig_l50 = as.numeric(assem_data$assembly_stats$contig_l50),
        contig_n50 = as.numeric(assem_data$assembly_stats$contig_n50),
        coverage = as.numeric(sub(assem_data$assembly_stats$genome_coverage, pattern = 'x$', replacement = '')),
        number_of_component_sequences = as.numeric(assem_data$assembly_stats$number_of_component_sequences),
        number_of_contigs = as.numeric(assem_data$assembly_stats$number_of_contigs),
        total_ungapped_length = as.numeric(assem_data$assembly_stats$total_ungapped_length),
        total_sequence_length = as.numeric(assem_data$assembly_stats$total_sequence_length),
        source_database = assem_data$source_database,
        is_type = "type_material" %in% names(assem_data),
        is_annotated = "annotation_info" %in% names(assem_data),
        is_atypical = "atypical" %in% names(assem_data$assembly_info),
        checkm_completeness = assem_data$checkm_info$completeness,
        checkm_contamination = assem_data$checkm_info$contamination
    )
    data_parts[sapply(data_parts, length) == 0 | sapply(data_parts, is.null)] <- NA
    output <- as.data.frame(data_parts)
    write.table(output, file = out_handle, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
}
close(out_handle)
close(in_handle)
