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

if (file.info(in_path)$size == 0) {
    # Handle empty files with only header
    write.table(header, file = out_path, sep = '\t', quote = FALSE, row.names = FALSE)
    quit(status = 0)
}
reports <- RcppSimdJson::fload(in_path, query = "/reports", always_list = TRUE, max_simplify_lvl = "list")[[1]]

# Handle null fields
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Get value helper
get_field <- function(x, path, default = NA_character_) {
  for (p in path) {
    if (is.null(x)) return(default)
    x <- x[[p]]
  }
  x %||% default
}

# Get host list
get_hosts <- function(r) {
  attrs <- r$assembly_info$biosample$attributes
  if (is.null(attrs)) return(NA_character_)
  hosts <- paste0(attrs$value[attrs$name == 'host'], collapse = ';')
  if (hosts == '') NA_character_ else hosts
}

num <- function(r, path) as.numeric(get_field(r, path, NA_character_))

out <- data.frame(
  accession       = vapply(reports, function(r) r$accession, character(1)),
  assembly_level  = vapply(reports, get_field, character(1), path = c("assembly_info", "assembly_level")),
  assembly_status = vapply(reports, get_field, character(1), path = c("assembly_info", "assembly_status")),
  assembly_type   = vapply(reports, get_field, character(1), path = c("assembly_info", "assembly_type")),
  hosts           = vapply(reports, get_hosts, character(1)),
  organism_name   = vapply(reports, function(r) get_field(r, c("organism", "organism_name")), character(1)),
  tax_id          = vapply(reports, function(r) as.character(get_field(r, c("organism", "tax_id"))), character(1)),
  contig_l50      = vapply(reports, num, double(1), path = c("assembly_stats", "contig_l50")),
  contig_n50      = vapply(reports, num, double(1), path = c("assembly_stats", "contig_n50")),
  coverage        = vapply(reports, function(r) {
                      v <- get_field(r, c("assembly_stats", "genome_coverage"))
                      if (is.na(v)) {
                          NA_real_
                      } else {
                          as.numeric(sub('x$', '', v))
                      }
                    }, double(1)),
  number_of_component_sequences = vapply(reports, num, double(1), path = c("assembly_stats", "number_of_component_sequences")),
  number_of_contigs             = vapply(reports, num, double(1), path = c("assembly_stats", "number_of_contigs")),
  total_ungapped_length         = vapply(reports, num, double(1), path = c("assembly_stats", "total_ungapped_length")),
  total_sequence_length         = vapply(reports, num, double(1), path = c("assembly_stats", "total_sequence_length")),
  source_database = vapply(reports, function(r) get_field(r, "source_database"), character(1)),
  is_type         = vapply(reports, function(r) "type_material" %in% names(r), logical(1)),
  is_annotated    = vapply(reports, function(r) "annotation_info" %in% names(r), logical(1)),
  is_atypical     = vapply(reports, function(r) "atypical" %in% names(r$assembly_info), logical(1)),
  checkm_completeness  = vapply(reports, num, double(1), path = c("checkm_info", "completeness")),
  checkm_contamination = vapply(reports, num, double(1), path = c("checkm_info", "contamination")),
  stringsAsFactors = FALSE
)

out$reference_id <- gsub('[\\/:*?"<>| .]', '_', out$accession)
out <- out[, c("reference_id", setdiff(names(out), "reference_id"))]
out$organism_name <- gsub('\\[|\\]', '', out$organism_name)

write.table(out, file = out_path, append = TRUE, sep = '\t', col.names = FALSE, row.names = FALSE, quote = FALSE)
