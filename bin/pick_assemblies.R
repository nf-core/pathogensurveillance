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


# Parse taxonomy inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     "~/projects/pathogensurveillance/work/5b/29ffb0e1e1add90d02de99307cffea/LF1_taxa_found.tsv",
#     "1",
#     "3",
#     "3",
#     "false",
#     "deleteme",
#     list.files("~/projects/pathogensurveillance/work/5b/29ffb0e1e1add90d02de99307cffea", pattern = '^[0-9]+.tsv$', full.names = TRUE)
# )

args <- as.list(args)
taxa_found_data <- read.table(args[[1]], header = TRUE, sep = '\t', comment.char = '')
n_ref_strains <- args[[2]]
n_ref_species <- args[[3]]
n_ref_genera <- args[[4]]
only_binomial <- as.logical(args[[5]])
out_name <- args[[6]]

# Parse input TSVs
if (length(args) < 7) {
    stop('No family-level reference metadata files supplied. Check input data.')
}
tsv_paths <- unlist(args[7:length(args)])
assem_data <- do.call(rbind, lapply(tsv_paths, function(path) {
    out <- read.table(path, header = TRUE, sep = '\t', comment.char = '')
    family_id <- gsub(basename(path), pattern = '.tsv', replacement = '', fixed = TRUE)
    out$family <- taxa_found_data$name[taxa_found_data$taxon_id == family_id]
    return(out)
}))

# Add column for modified ID
modified_id <- gsub(assem_data$accession, pattern = '[\\/:*?"<>| .]', replacement = '_')
assem_data <- cbind(reference_id = modified_id, assem_data)

# Add taxon info columns
assem_data$organism_name <- gsub(assem_data$organism_name, pattern = '[', replacement = '', fixed = TRUE)
assem_data$organism_name <- gsub(assem_data$organism_name, pattern = ']', replacement = '', fixed = TRUE)
genus_prefixes <- c('Candidatus')
genus_prefix_pattern <- paste0(genus_prefixes, collapse = '|')
binomial_pattern <- paste0('^((?:', genus_prefix_pattern, ')?) ?([a-zA-Z0-9.]+) ?([a-zA-Z0-9.]+) ?(.*)$')
assem_data$species <- trimws(gsub(assem_data$organism_name, pattern = binomial_pattern, replacement = '\\1 \\2 \\3', ignore.case = TRUE))
assem_data$genus <- trimws(gsub(assem_data$organism_name, pattern = binomial_pattern, replacement = '\\1 \\2', ignore.case = TRUE))

# Filter out references with non-standard names
is_ambiguous <- function(x) {
    ambiguous_words <- c(
        'uncultured',
        'unknown',
        'incertae sedis',
        'sp\\.',
        'cf\\.',
        'endosymbiont',
        'symbiont',
        'bacterium'
    )
    ambiguous_pattern <- paste0('\\b', ambiguous_words, '\\b', collapse = '|')
    grepl(x, pattern = ambiguous_pattern, ignore.case = TRUE)
}
is_latin_binomial <- function(x) {
    grepl(x, pattern = '^[a-zA-Z]+ [a-zA-Z]+($| ).*$') & ! is_ambiguous(x)
}
if (only_binomial) {
    assem_data <- assem_data[is_latin_binomial(assem_data$species), ]
}

# Parse "count" arguments which can be a number or a percentage
get_count <- function(choices, count) {
    if (grepl(count, pattern = "%$")) {
        prop <- as.numeric(sub(count, pattern = "%", replacement = "")) / 100
        count <- ceiling(choices * prop)
        return(min(c(choices, count)))
    } else {
        count <- as.numeric(count)
        return(min(c(choices, count)))
    }
}

# Sort references by desirability
priority <- order(decreasing = TRUE,
    assem_data$is_atypical == FALSE,
    assem_data$is_type, # Is type strain
    assem_data$source_database == 'SOURCE_DATABASE_REFSEQ', # Is a RefSeq reference
    is_latin_binomial(assem_data$species), # Has a species epithet
    assem_data$is_annotated,
    factor(assem_data$assembly_level, levels = c("Contig", "Scaffold", "Chromosome", "Complete Genome"), ordered = TRUE),
    assem_data$checkm_completeness,
    -1 * assem_data$checkm_contamination,
    assem_data$contig_l50,
    assem_data$coverage
)
assem_data <- assem_data[priority, ]

# Initialize column to hold which level an assembly is selected for
assem_data$selection_rank <- NA
assem_data$selection_taxon <- NA
assem_data$selection_subtaxon <- NA

# Select representatives for each rank
select_for_rank <- function(assem_data, query_taxa, rank, subrank, count_per_rank, count_per_subrank = 1)  {
    for (tax in query_taxa) {
        # Get assembly indexes for every subtaxon
        tax_to_consider <- (assem_data[[rank]] == tax | assem_data[[subrank]] == tax) & is.na(assem_data$selection_rank)
        subtaxa_found <- unique(assem_data[[subrank]][tax_to_consider])
        subtaxa_found <- subtaxa_found[! subtaxa_found %in% c(assem_data$selection_taxon, assem_data$selection_subtaxon)] # Dont include the data for taxa already chosen
        selected <- lapply(subtaxa_found, function(subtax) {
            which(assem_data[[subrank]] == subtax & is.na(assem_data$selection_rank))
        })
        names(selected) <- subtaxa_found

        # Parse count attributes, which can be percentages or integers
        count_per_rank <- get_count(length(selected), count_per_rank)
        count_per_subrank <- get_count(length(selected), count_per_subrank)

        # Pick subtaxa with the most assemblies and best mean attributes (based on order in input)
        mean_index <- vapply(selected, mean, FUN.VALUE = numeric(1))
        subtaxa_count <- vapply(selected, length, FUN.VALUE = numeric(1))
        selection_priority <- order(decreasing = TRUE,
            is_ambiguous(names(selected)) == FALSE,
            subtaxa_count,
            -mean_index
        )
        selected <- selected[selection_priority]
        selected <- selected[seq_len(min(c(count_per_rank, length(selected))))]

        # Pick representatives of subtaxa with best attributes (based on order in input)
        selected <- lapply(selected, function(x) {
            x[seq_len(min(c(count_per_subrank, length(x))))]
        })

        # Record data on selected assemblies
        selected <- unlist(selected)
        assem_data$selection_rank[selected] <- rank
        assem_data$selection_taxon[selected] <- tax
        assem_data$selection_subtaxon[selected] <- assem_data[[subrank]][selected]
    }
    return(assem_data)
}
assem_data <- select_for_rank(
    assem_data,
    query_taxa = taxa_found_data$name[taxa_found_data$rank == 'species'],
    rank = 'species',
    subrank = 'organism_name',
    count_per_rank = n_ref_strains
)
assem_data <- select_for_rank(
    assem_data,
    query_taxa = taxa_found_data$name[taxa_found_data$rank == 'genus'],
    rank = 'genus',
    subrank = 'species',
    count_per_rank = n_ref_species
)
assem_data <- select_for_rank(
    assem_data,
    query_taxa = taxa_found_data$name[taxa_found_data$rank == 'family'],
    rank = 'family',
    subrank = 'genus',
    count_per_rank = n_ref_genera
)

result <- assem_data[! is.na(assem_data$selection_taxon), ]

# Reformat results to the same format as the user-defined metadata
if (nrow(result) == 0) {
    formatted_result <- data.frame(
        ref_id = character(0),
        ref_name = character(0),
        ref_description = character(0),
        ref_path = character(0),
        ref_ncbi_accession = character(0),
        ref_ncbi_query = character(0),
        ref_ncbi_query_max = character(0),
        ref_primary_usage = character(0),
        ref_contextual_usage = character(0),
        ref_color_by = character(0)
    )
} else {
    suffix <- paste0(
        ifelse(result$is_type, 'T', ''),
        ifelse(result$source_database == 'SOURCE_DATABASE_REFSEQ', 'R', ''),
        ifelse(result$is_atypical, 'A', '')
    )
    formatted_result <- data.frame(
        ref_id = result$reference_id,
        ref_name = result$organism_name,
        ref_description = paste0(
            result$organism_name, ' ',
            result$accession,
            ifelse(nchar(suffix) > 0, paste0(' ', suffix), '')
        ),
        ref_path = '',
        ref_ncbi_accession = result$accession,
        ref_ncbi_query = '',
        ref_ncbi_query_max = '',
        ref_primary_usage = 'optional',
        ref_contextual_usage = 'optional',
        ref_color_by = ''
    )
}

# Save to output file
write.table(result, file = paste0(out_name, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(formatted_result, file = paste0(out_name, '_formatted.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
