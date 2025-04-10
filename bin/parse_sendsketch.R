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


# Options
ani_threshold <- c(species = 90, genus = 80, family = 50)  # These numbers are total guesses. TODO: find reasonable defaults (issue #11)
complt_threshold <- c(species = 40, genus = 15, family = 5) # These numbers are total guesses. TODO: find reasonable defaults (issue #11)

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- list('/home/fosterz/projects/pathogensurveillance/work/12/25ffa8ac8801703ed6717d3f3f43fa/LF1.txt')
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

