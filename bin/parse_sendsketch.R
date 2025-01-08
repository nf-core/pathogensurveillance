#!/usr/bin/env Rscript

# Options
ani_threshold <- c(species = 90, genus = 80, family = 50)  # These numbers are total guesses. TODO: find reasonable defaults (issue #11)
complt_threshold <- c(species = 40, genus = 15, family = 5) # These numbers are total guesses. TODO: find reasonable defaults (issue #11)

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- list('~/projects/pathogensurveillance/work/66/c64114197ddf5ea301f6cd3ffc4e23/SRR26197583.txt')
data <- read.csv(args[[1]], skip = 2, header = TRUE, sep = '\t')

# Format table
data$ANI <- as.numeric(sub(pattern = "%", replacement = "", fixed = TRUE, data$ANI))
data$Complt <- as.numeric(sub(pattern = "%", replacement = "", fixed = TRUE, data$Complt))

# Filter data by threshold and extract passing taxon names
filter_and_extract <- function(table, level, rank_code) {
    # Filter table by ANI threshold and completeness
    table <- table[table$ANI > ani_threshold[level] & table$Complt > complt_threshold[level], ]
    # If no taxa are found, return an empty vector
    if (nrow(table) == 0) {
        return(character())
    }
    # Make list of vectors of taxa named by rank
    parsed_classifications <- lapply(strsplit(table$taxonomy, split = ';', fixed = TRUE), function(split_text) {
        rank_and_taxon <- strsplit(split_text, split = ':', fixed = TRUE)
        output <- unlist(lapply(rank_and_taxon, `[`, 2))
        names(output) <- unlist(lapply(rank_and_taxon, `[`, 1))
        return(output)
    })
    # Extract taxa for rank of interest
    output <- unname(unlist(lapply(parsed_classifications, `[`, rank_code)))
    # Remove NAs and make unique
    output <- unique(output[!is.na(output)])
    return(output)
}

# Extract taxon info
genus <- filter_and_extract(data, "genus", "g")
family <- filter_and_extract(data, "family", "f")
species <- filter_and_extract(data, "species", "s")

# Write output
writeLines(species, "species.txt")
writeLines(genus, "genera.txt")
writeLines(family, "families.txt")

# Save kingdom
kingdom <- gsub("sk:(.+?);[kp]:.*", "\\1", data$taxonomy[1])
writeLines(kingdom, "kingdom.txt")

# Save full taxonomic classification
writeLines(data$taxonomy[1], "classification.txt")

