#!/usr/bin/env Rscript

# Options
ani_threshold <- c(species = 90, genus = 90, family = 70)  # These numbers are total guesses. TODO: find reasonable defaults (issue #11)
complt_threshold <- c(species = 60, genus = 40, family = 30) # These numbers are total guesses. TODO: find reasonable defaults (issue #11)

# Parse inputs
args = commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], skip = 2, header = TRUE, sep = '\t')

# Format table
data$ANI <- as.numeric(sub(pattern = "%", replacement = "", fixed = TRUE, data$ANI))
data$Complt <- as.numeric(sub(pattern = "%", replacement = "", fixed = TRUE, data$Complt))

# Filter data by threshold and extract passing taxon names
filter_and_extract <- function(table, level, pattern) {
  table <- table[table$ANI > ani_threshold[level] & table$Complt > complt_threshold[level], ]
  gsub(pattern, "\\1", table$taxonomy)
}

# Extract taxon info
genus <- filter_and_extract(data, "genus", ".*g:(.+);s:.*")
family <- filter_and_extract(data, "family", ".*f:(.+);g:.*")
species <- filter_and_extract(data, "species", ".*;s:([a-zA-Z0-9 .]+);?.*")

# Write output
writeLines(unique(species), "species.txt")
writeLines(unique(genus), "genera.txt")
writeLines(unique(family), "families.txt")

# Save kingdom 
kingdom <- gsub("sk:(.+);p:.*", "\\1", data$taxonomy[1])
writeLines(kingdom, "kingdom.txt")

# Save full taxonomic classification
writeLines(data$taxonomy[1], "classification.txt")

