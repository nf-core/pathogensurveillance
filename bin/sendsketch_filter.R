#!/usr/bin/env -S Rscript --vanilla

#args=c("C:\Users\smile\Downloads\work\22-299_T1.txt")

# Options
ani_threshold <- 90

# Parse inputs
args = commandArgs(trailingOnly = TRUE)
data <- read.csv(args[1], skip = 2, header = TRUE, sep = '\t')

# Format table
data$ANI <- as.numeric(sub(pattern = "%", replacement = "", fixed = TRUE, data$ANI))

# Filter table by match quality
data <- data[data$ANI > ani_threshold, ]

# Extract taxon info
data$genus <- gsub(".*g:(.+);s:.*", "\\1", data$taxonomy)
data$family <- gsub(".*f:(.+);g:.*", "\\1", data$taxonomy)
data$species <- gsub(".*;s:([a-zA-Z0-9 .]+);?.*", "\\1", data$taxonomy)

# Write output
writeLines(unique(data$species), "species.txt")
writeLines(unique(data$genus), "genera.txt")
writeLines(unique(data$family), "families.txt")

# Save kingdom 
kingdom <- gsub("sk:(.+);p:.*", "\\1", data$taxonomy[1])
writeLines(kingdom, "kingdom.txt")

# Save full taxonomic classification
writeLines(data$taxonomy[1], "classification.txt")
