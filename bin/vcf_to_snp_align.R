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

# Constants
missing_data_char <- '-'
max_missing_data_prop <- 0.5  # SNPs are removed if the proportion of samples with data is greater than this. Set to 1 or NA to disable SNP filtering.
remove_missing_samples <- TRUE # Whether to remove samples with only missing data
pick_best_allele_for_haploids <- FALSE # Use depth information to pick the best allele for haploids that have information for two "alleles". Slower.
missing_sample_file_path <- 'removed_sample_ids.txt'

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#    '~/projects/pathogensurveillance/work/09/440702826216a7101811d85a080794/subgroup_GCF_903978215_1.vcffilter.vcf.gz',
#    '~/projects/pathogensurveillance/work/09/440702826216a7101811d85a080794/subgroup_GCF_903978215_1.tsv',
#    '1000000',
#    'deleteme.fasta'
# )
names(args) <- c('vcf_path', 'ploidy_data_path', 'max_variants', 'out_path')
args <- as.list(args)

# Get header from VCF file
file_handle <- file(args$vcf_path, "r")
while ( TRUE ) {
    line <- readLines(file_handle, n = 1)
    if ( length(line) == 0 ) {
        stop('Cannot find VCF header.')
    }
    if (grepl(line, pattern = '^#CHROM')) {
        header_line <- line
        break
    }
}
no_variants <- length(readLines(file_handle, n = 1)) == 0
close(file_handle)
header <- strsplit(header_line, split = '\t')[[1]]
header[1] <- 'CHROM'

# Read VCF file and name by header
if (no_variants) {
    vcf_data <- as.data.frame(rep(list(character(0)), length(header)))
} else {
    vcf_data <- read.delim(file = args$vcf_path, sep = '\t', comment.char = '#', header = FALSE, nrows = as.numeric(args$max_variants))
}
colnames(vcf_data) <- header

# Read ploidy data file
ploidy_data <- read.csv(args$ploidy_data_path, sep = '\t')
ploidy_key <- stats::setNames(ploidy_data$ploidy, ploidy_data$mapping_id)

# Make table with just alleles for each sample
sample_ids <- colnames(vcf_data)[10:length(colnames(vcf_data))]
gt_index <- vapply(strsplit(vcf_data$FORMAT, split = ':'), FUN.VALUE = integer(1), function(x) {
    which(x == 'GT')
})
ad_index <- vapply(strsplit(vcf_data$FORMAT, split = ':'), FUN.VALUE = integer(1), function(x) {
    which(x == 'AD')
})
vcf_data[sample_ids] <- lapply(sample_ids, function(samp_id) {
    raw_values <- vcf_data[[samp_id]]
    split_values <- strsplit(raw_values, split = ':')
    genotypes <- mapply(`[`, split_values, gt_index, SIMPLIFY = TRUE)
    if (ploidy_key[samp_id] == 1 && pick_best_allele_for_haploids) {
        genotypes <- strsplit(genotypes, '/')
        depths <- vapply(seq_along(split_values), FUN.VALUE = character(1), function(i) split_values[[i]][ad_index[i]])
        depths <- strsplit(depths, ',')
        depths <- lapply(depths, function(x) as.numeric(x))
        depths <- lapply(seq_along(depths), function(i) {
            genos <- genotypes[[i]]
            genos[genos == '.'] <- NA
            depths[[i]][as.numeric(genos) + 1] # Is this + 1 supposed to be here?
        })
        is_all_na <- vapply(depths, FUN.VALUE = logical(1), function(x) all(is.na(x)))
        genotypes <- vapply(seq_along(genotypes), FUN.VALUE = character(1), function(i) {
            if (all(is.na(depths[[i]]))) {
                return('.')
            } else {
                return(as.character(genotypes[[i]][which.max(depths[[i]])]))
            }
        })
    }
    return(genotypes)
})


# Create IUPAC key with all combinations precomputed for speed
iupac_key <- c(
    A='A',
    C='C',
    G='G',
    T='T',
    AC='M',
    AG='R',
    AT='W',
    CG='S',
    CT='Y',
    GT='K',
    ACG='V',
    ACT='H',
    AGT='D',
    CGT='B',
    ACGT='N'
)
permute_char <- function(char) {
    if (nchar(char) == 1) {
        return(char)
    }
    split <- strsplit(char, split = '')[[1]]
    options <- expand.grid(rep(list(split), length(split)))
    is_duplicated <- apply(options, MARGIN = 1, function(x) any(duplicated(x)))
    unname(apply(options[! is_duplicated, ], MARGIN = 1, paste0, collapse = ''))
}
iupac_permutations <- lapply(names(iupac_key), permute_char)
expanded_iupac_key <- rep(iupac_key, vapply(iupac_permutations, FUN.VALUE = numeric(1), length))
names(expanded_iupac_key) <- unlist(iupac_permutations)

# Convert allele numeric codes to AGCT + IUPAC ambiguity codes
vcf_data[sample_ids] <- lapply(sample_ids, function(samp_id) {
    if (ploidy_key[samp_id] == 1 && ! pick_best_allele_for_haploids) {
        # Parse indexes of alleles for each haplotype
        allele_codes <- substring(vcf_data[[samp_id]], 1, 1) # A hack for speed. Ignores alternative allele and will return incorrect results if there are more than 9 alleles
        allele_codes[allele_codes == '.'] <- NA
        allele_codes <- as.numeric(allele_codes)

        # Convert allele codes to sequence
        allele_code_key <- mapply(c, vcf_data$REF, strsplit(vcf_data$ALT, split = ','), SIMPLIFY = FALSE)
        alleles <- mapply(allele_codes + 1, allele_code_key, SIMPLIFY = TRUE, FUN = function(a, k) {
            expanded_iupac_key[k[a]]
        })
        return(unname(alleles))
    } else {
        # Parse indexes of alleles for each haplotype
        allele_codes <- strsplit(vcf_data[[samp_id]], split = '/')
        allele_codes <- lapply(allele_codes, function(x) as.numeric(x[x != '.']))

        # Convert allele codes to sequence
        allele_code_key <- mapply(c, vcf_data$REF, strsplit(vcf_data$ALT, split = ','), SIMPLIFY = FALSE)
        alleles <- mapply(allele_codes, allele_code_key, SIMPLIFY = FALSE, FUN = function(a, k) {
            out <- k[a + 1]
            if (any(nchar(out) != 1)) { # Return NA if not a simple SNP
                return(NA_character_)
            } else {
                return(expanded_iupac_key[paste0(unique(out), collapse = '')]) # Convert alleles to a single sequence, using IUPAC ambiguity codes as needed
            }
        })
        return(unname(unlist(alleles)))
    }
})

# Remove samples with only missing data
if (remove_missing_samples) {
    is_only_missing_data <- vapply(vcf_data[sample_ids], FUN.VALUE = logical(1), function(x) all(is.na(x)))
    samples_to_remove <- sample_ids[is_only_missing_data]
    vcf_data[samples_to_remove] <- NULL
    if (all(is_only_missing_data)) {
        vcf_data <- vcf_data[numeric(0), , drop = FALSE]
    }
    sample_ids <- sample_ids[! sample_ids %in% samples_to_remove]
    writeLines(samples_to_remove, missing_sample_file_path)
}

# Remove variants with missing data
if (! is.na(max_missing_data_prop) && max_missing_data_prop < 1) {
    prop_missing <- apply(vcf_data[sample_ids], MARGIN = 1, simplify = TRUE, function(x) sum(is.na(x)) / length(x))
    vcf_data <- vcf_data[prop_missing <= max_missing_data_prop, , drop = FALSE]
}

# Convert missing data to a character
vcf_data[sample_ids] <- lapply(vcf_data[sample_ids], function(x) ifelse(is.na(x), missing_data_char, x))

# Remove indels relative to the reference
vcf_data <- vcf_data[nchar(vcf_data$REF) == 1, , drop = FALSE]

# Convert to sequences
reference_seq <- paste0(vcf_data$REF, collapse = '')
names(reference_seq) <- "REF"
sample_seqs <- vapply(vcf_data[sample_ids], FUN.VALUE = character(1), paste0, collapse = '')
seqs <- c(reference_seq, sample_seqs)

# Write output FASTA file
writeLines(paste0('>', names(seqs), '\n', seqs), args$out_path)
