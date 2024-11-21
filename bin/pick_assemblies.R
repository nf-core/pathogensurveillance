#!/usr/bin/env Rscript

library(jsonlite)

# Set random number generator seed
set.seed(1)

# Parse taxonomy inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/SRR29399310_families.txt',
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/SRR29399310_genera.txt',
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/SRR29399310_species.txt',
#     '5', '10', '10', 'false', 'output.tsv',
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/Saprolegniaceae.json',
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/Peronosporaceae.json',
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/Pythiaceae.json',
#     '/home/fosterz/projects/pathogensurveillance/work/cc/3bffaedde209f7e754c9c9dc960994/Staphylococcaceae.json'
# )
args <- as.list(args)
families <- readLines(args[[1]])
genera <- readLines(args[[2]])
species <- readLines(args[[3]])
n_ref_strains <- args[[4]]
n_ref_species <- args[[5]]
n_ref_genera <- args[[6]]
only_binomial <- as.logical(args[[7]])
out_path <- args[[8]]

# Parse input JSONs
if (length(args) < 9) {
    stop('No family-level reference metadata files supplied. Check input data.')
}
json_paths <- unlist(args[9:length(args)])
json_families <- gsub(basename(json_paths), pattern = '.json', replacement = '', fixed = TRUE)
json_data <- lapply(seq_along(json_paths), function(index) {
    parsed_json <- lapply(readLines(json_paths[index]), jsonlite::fromJSON)
    output <- do.call(rbind, lapply(parsed_json, function(assem_data) {
        attributes <- assem_data$assembly_info$biosample$attributes
        hosts <- paste0(attributes$value[attributes$name == 'host'], collapse = ';')
        data.frame(
            accession = assem_data$accession,
            assembly_level = assem_data$assembly_info$assembly_level,
            assembly_status = assem_data$assembly_info$assembly_status,
            assembly_type = assem_data$assembly_info$assembly_type,
            hosts = ifelse(hosts == '', NA_character_, hosts),
            organism_name = assem_data$organism$organism_name,
            tax_id = as.character(assem_data$organism$tax_id),
            contig_l50 = as.numeric(assem_data$assembly_stats$contig_l50),
            contig_n50 = as.numeric(assem_data$assembly_stats$contig_n50),
            number_of_component_sequences = as.numeric(assem_data$assembly_stats$number_of_component_sequences),
            number_of_contigs = as.numeric(assem_data$assembly_stats$number_of_contigs),
            total_ungapped_length = as.numeric(assem_data$assembly_stats$total_ungapped_length),
            total_sequence_length = as.numeric(assem_data$assembly_stats$total_sequence_length),
            source_database = assem_data$source_database,
            is_type = "type_material" %in% names(assem_data),
            is_annotated = "annotation_info" %in% names(assem_data),
            is_atypical = "atypical" %in% names(assem_data$assembly_info)
        )
    }))
    if (!is.null(output)) {
        output$family <- rep(json_families[index], nrow(output))
    }
    return(output)
})
assem_data <- do.call(rbind, json_data)
if (is.null(assem_data)) {
    assem_data <- data.frame(
        accession = character(0),
        assembly_level = character(0),
        assembly_status = character(0),
        assembly_type = character(0),
        hosts = character(0),
        organism_name = character(0),
        tax_id = character(0),
        contig_l50 = numeric(0),
        contig_n50 = numeric(0),
        number_of_component_sequences = numeric(0),
        number_of_contigs = numeric(0),
        total_ungapped_length = numeric(0),
        total_sequence_length = numeric(0),
        source_database = character(0),
        is_type = logical(0),
        is_annotated = logical(0),
        is_atypical = logical(0)
    )
}

# Add column for modified ID
modified_id <- gsub(assem_data$accession, pattern = '[\\/:*?"<>| .]', replacement = '_')
assem_data <- cbind(reference_id = modified_id, assem_data)

# Add taxon info columns
assem_data$organism_name <- gsub(assem_data$organism_name, pattern = '[', replacement = '', fixed = TRUE)
assem_data$organism_name <- gsub(assem_data$organism_name, pattern = ']', replacement = '', fixed = TRUE)
assem_data$species <- gsub(assem_data$organism_name, pattern = '([a-zA-Z0-9.]+) ([a-zA-Z0-9.]+) (.*)', replacement = '\\1 \\2')
assem_data$genus <- gsub(assem_data$organism_name, pattern = '([a-zA-Z0-9.]+) (.*)', replacement = '\\1')

# Filter out references with non-standard names
is_latin_binomial <- function(x) {
    grepl(x, pattern = '^[a-zA-Z]+ [a-zA-Z]+($| ).*$')
}
if (only_binomial) {
    assem_data <- assem_data[is_latin_binomial(assem_data$species), ]
}

# Parse "count" arguments which can be a number or a percentage
get_count <- function(table, count) {
    if (grepl(count, pattern = "%$")) {
        prop <- as.numeric(sub(count, pattern = "%", replacement = "")) / 100
        count <- ceiling(nrow(table) * prop)
        return(min(c(nrow(table), count)))
    } else {
        count <- as.numeric(count)
        return(min(c(nrow(table), count)))
    }
}

# Function to select best references given metadata.
#  This is used to pick between multiple references that meet a taxonomic criteria.
#  The selection is based on assembly quality and the species name
select_references <- function(ref_data, n = 1) {
    priority <- order(
        decreasing = TRUE,
        ref_data$is_atypical == FALSE,
        ref_data$is_type, # Is type strain
        ref_data$source_database == 'SOURCE_DATABASE_REFSEQ', # Is a RefSeq reference
        is_latin_binomial(ref_data$species), # Has a species epithet
        ref_data$is_annotated,
        factor(ref_data$assembly_level, levels = c("Contig", "Scaffold", "Chromosome", "Complete Genome"), ordered = TRUE), 
        ref_data$contig_n50
    )
    return(ref_data[priority[seq_len(n)], ])
}

# Pick representatives for each species
sp_stats <- assem_data[assem_data$species %in% species, , drop = FALSE]
sp_stats <- lapply(split(sp_stats, sp_stats$species), function(per_sp_data) {
    supset_per_sp_data <- lapply(split(per_sp_data, per_sp_data$organism_name), function(per_org_data) {
        select_references(per_org_data, n = 1)
    })
    supset_per_sp_data <- do.call(rbind, supset_per_sp_data)
    select_references(supset_per_sp_data, n = get_count(supset_per_sp_data, n_ref_strains))
})
if (length(sp_stats) == 0 )  {
    sp_stats <- assem_data[numeric(0), ]
} else {
    sp_stats <- do.call(rbind, sp_stats)
    rownames(sp_stats) <- NULL
} 

# Pick representatives for each genus
gn_stats <- assem_data[assem_data$genus %in% genera, , drop = FALSE]
gn_stats <- lapply(split(gn_stats, gn_stats$genus), function(per_gn_data) {
    subset_per_gn_data <- lapply(split(per_gn_data, per_gn_data$species), function(per_org_data) {
        select_references(per_org_data, n = 1)
    })
    subset_per_gn_data <- do.call(rbind, subset_per_gn_data)
    subset_per_gn_data <- subset_per_gn_data[! subset_per_gn_data$species %in% sp_stats$species, ] # Dont include the species already chosen
    select_references(subset_per_gn_data, n = get_count(subset_per_gn_data, n_ref_species))
})
if (length(gn_stats) == 0 ) {
    gn_stats <- assem_data[numeric(0), ]
} else {
    gn_stats <- do.call(rbind, gn_stats)
    rownames(gn_stats) <- NULL
}

# Pick representatives for each family
fa_stats <- assem_data[assem_data$family %in% families, , drop = FALSE]
fa_stats <- lapply(split(fa_stats, fa_stats$family), function(per_fm_data) {
    subset_per_fm_data <- lapply(split(per_fm_data, per_fm_data$genus), function(per_gn_data) {
        select_references(per_gn_data, n = 1)
    })
    subset_per_fm_data <- do.call(rbind, subset_per_fm_data)
    subset_per_fm_data <- subset_per_fm_data[! subset_per_fm_data$genus %in% gn_stats$genus, ] # Dont include the genera already chosen
    select_references(subset_per_fm_data, n = get_count(subset_per_fm_data, n_ref_genera))
})
if (length(fa_stats) == 0 ) {
    fa_stats <- assem_data[numeric(0), ]
} else {
    fa_stats <- do.call(rbind, fa_stats)
    rownames(fa_stats) <- NULL
}

# Combine results
result <- rbind(sp_stats, gn_stats, fa_stats)

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
    formatted_result <- data.frame(
        ref_id = result$reference_id,
        ref_name = result$organism_name,
        ref_description = paste0(result$organism_name, ' (', 
                                 result$accession,
                                 ifelse(result$is_type, '; Type', ''),
                                 ifelse(result$source_database == 'SOURCE_DATABASE_REFSEQ', '; RefSeq', ''),
                                 ifelse(result$is_atypical, '; Atypical', ''),
                                 ifelse(is.na(result$hosts), '', paste0('; Host: ', result$hosts)),
                                 ')'),
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
write.table(formatted_result, file = out_path, sep = '\t', quote = FALSE, row.names = FALSE)
write.table(assem_data, file = 'merged_assembly_stats.tsv', sep = '\t', quote = FALSE, row.names = FALSE)
