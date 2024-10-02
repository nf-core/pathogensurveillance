#!/usr/bin/env Rscript

# Set random number generator seed
set.seed(1)

# Options
min_coverage <- 30

# Parse taxonomy inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#      '/home/fosterz/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/ERR3500120_families.txt',
#      '/home/fosterz/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/ERR3500120_genera.txt',
#      '/home/fosterz/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/ERR3500120_species.txt',
#      '30', '20', '10',
#      'ERR3500120.tsv',
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Aleyrodidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Amborellaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Aphididae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Buthidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Cerambycidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Cercopithecidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Cheirogaleidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Chrysomelidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Cordycipitaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Cricetidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Cucurbitaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Culicidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Dasypodidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Dasyuridae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Delphacidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Fabaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Formicidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Galagidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Halomonadaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Lampyridae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Liviidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Malvaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Micrococcaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Muscidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Nectriaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Nitidulidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Orchidaceae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Otariidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Pentatomidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Tephritidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Theridiidae.tsv', 
#      '~/files/projects/pathogensurveillance/work/98/9b5787f3c655cd5f51f2194df064cf/Xanthomonadaceae.tsv' )
args <- as.list(args)
families <- readLines(args[[1]])
genera <- readLines(args[[2]])
species <- readLines(args[[3]])
n_ref_strains <- args[[4]]
n_ref_species <- args[[5]]
n_ref_genera <- args[[6]]
out_path <- args[[7]]

# Parse input TSVs
if (length(args) < 8) {
    stop('No family-level reference metadata files supplied. Check input data.')
}
tsv_paths <- unlist(args[8:length(args)])
tsv_families <- gsub(basename(tsv_paths), pattern = '.tsv', replacement = '', fixed = TRUE)
tsv_data <- lapply(seq_along(tsv_paths), function(index) {
    output <- read.csv(tsv_paths[index], sep = '\t')
    output$family <- rep(tsv_families[index], nrow(output))
    return(output)
})
assem_data <- do.call(rbind, tsv_data)
assem_data$Coverage <- as.numeric(assem_data$Coverage)
assem_data$ScaffoldN50 <- as.numeric(assem_data$ScaffoldN50)
assem_data$genus <- gsub(assem_data$SpeciesName, pattern = '([a-zA-Z0-9.]+) (.*)', replacement = '\\1')

# Add column for modified ID
modified_id <- gsub(assem_data$LastMajorReleaseAccession, pattern = '[\\/:*?"<>| .]', replacement = '_')
assem_data <- cbind(reference_id = modified_id, assem_data)

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

# Quality control
assem_data <- assem_data[! is.na(assem_data$Coverage) & assem_data$Coverage >= min_coverage, , drop = FALSE]

# Pick representatives for each species
sp_stats <- assem_data[assem_data$SpeciesName %in% species, , drop = FALSE]
sp_stats <- lapply(split(sp_stats, sp_stats$SpeciesName), function(per_sp_data) {
  lapply(split(per_sp_data, per_sp_data$Organism), function(per_org_data) {
    per_org_data[which.max(per_org_data$ScaffoldN50), ]
  })
})
if (length(sp_stats) == 0 )  {
  sp_stats <- assem_data[numeric(0), ]
} else {
  sp_stats <- do.call(rbind, unlist(sp_stats, recursive = FALSE))
  rownames(sp_stats) <- NULL
  sp_stats <- sp_stats[order(sp_stats$ScaffoldN50, decreasing = TRUE)[seq_len(get_count(sp_stats, n_ref_strains))], ]
} 

# Pick representatives for each genus
gn_stats <- assem_data[assem_data$genus %in% genera, , drop = FALSE]
gn_stats <- lapply(split(gn_stats, gn_stats$genus), function(per_gn_data) {
  lapply(split(per_gn_data, per_gn_data$SpeciesName), function(per_org_data) {
    per_org_data[which.max(per_org_data$ScaffoldN50), ]
  })
})
if (length(gn_stats) == 0 ) {
  gn_stats <- assem_data[numeric(0), ]
} else {
  gn_stats <- do.call(rbind, unlist(gn_stats, recursive = FALSE))
  rownames(gn_stats) <- NULL
  gn_stats <- gn_stats[! gn_stats$SpeciesName %in% sp_stats$SpeciesName, ] # Dont include the species already chosen
  gn_stats <- gn_stats[order(gn_stats$ScaffoldN50, decreasing = TRUE)[seq_len(get_count(gn_stats, n_ref_species))], ]
}

# Pick representatives for each family
fa_stats <- assem_data[assem_data$family %in% families, , drop = FALSE]
fa_stats <- lapply(split(fa_stats, fa_stats$family), function(per_fm_data) {
  lapply(split(per_fm_data, per_fm_data$genus), function(per_gn_data) {
    per_gn_data[which.max(per_gn_data$ScaffoldN50), ]
  })
})
if (length(fa_stats) == 0 ) {
  fa_stats <- assem_data[numeric(0), ]
} else {
  fa_stats <- do.call(rbind, unlist(fa_stats, recursive = FALSE))
  rownames(fa_stats) <- NULL
  fa_stats <- fa_stats[! fa_stats$genus %in% gn_stats$genus, ] # Dont include the genera already chosen
  fa_stats <- fa_stats[order(fa_stats$ScaffoldN50, decreasing = TRUE)[seq_len(get_count(fa_stats, n_ref_genera))], ]
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
    ref_name = gsub(result$Organism, pattern = ' \\(.+\\)$', replacement = ''),
    ref_description = paste0(gsub(result$Organism, pattern = ' \\(.+\\)$', replacement = ''), ' (', result$LastMajorReleaseAccession, ')'),
    ref_path = '',
    ref_ncbi_accession = result$LastMajorReleaseAccession,
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
