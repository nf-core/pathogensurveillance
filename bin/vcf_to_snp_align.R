#!/usr/bin/env Rscript

# Constants
missing_data_char <- '-'
max_missing_data_prop <- 0.5  # SNPs are removed if the proportion of samples with data is greater than this. Set to 1 or NA to disable SNP filtering.
remove_missing_samples <- TRUE # Whether to remove samples with only missing data
missing_sample_file_path <- 'removed_sample_ids.txt'

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
# args <- c(
#    '~/projects/pathogensurveillance/work/7e/982f33311eb286918183ec32916767/mixed_GCF_042647405_1.vcffilter.vcf.gz',
#    '~/projects/pathogensurveillance/work/7e/982f33311eb286918183ec32916767/mixed_GCF_042647405_1.csv',
#    'deleteme.fasta'
# )
names(args) <- c('vcf_path', 'ploidy_data_path', 'out_path')
args <- as.list(args)

# Read VCF file
header_data <- readLines(args$vcf_path)
header_line <- header_data[grepl(header_data, pattern = '^#CHROM')]
header <- strsplit(header_line, split = '\t')[[1]]
header[1] <- 'CHROM'
vcf_data <- read.delim(file = args$vcf_path, sep = '\t', comment.char = '#', header = FALSE)
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
    vapply(seq_len(nrow(vcf_data)), FUN.VALUE = character(1), function(index) {
        genotypes <- strsplit(vcf_data[index, samp_id], split = ':')[[1]][gt_index[index]]
        genotypes <- strsplit(genotypes, '/')[[1]]
        genotypes[genotypes == '.'] <- NA
        genotypes <- as.numeric(genotypes)
        depths <- strsplit(vcf_data[index, samp_id], split = ':')[[1]][ad_index[index]]
        depths <- as.numeric(strsplit(depths, ',')[[1]])
        depths <- depths[genotypes]
        ploidy <- ploidy_key[samp_id]
        if (ploidy == 1) {
            if (all(is.na(depths))) {
                return('.')
            } else {
                return(as.character(genotypes[which.max(depths)]))
            }
        } else {
            genotypes <- ifelse(is.na(genotypes), '.', as.character(genotypes))
            return(paste0(genotypes, collapse = '/'))
        }
    })
})

# Convert allele numeric codes to AGCT + IUPAC ambiguity codes
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
vcf_data[sample_ids] <- lapply(sample_ids, function(samp_id) {
    vapply(seq_len(nrow(vcf_data)), FUN.VALUE = character(1), function(index) {
        # Parse indexes of alleles for each haplotype
        allele_codes <- strsplit(vcf_data[index, samp_id], split = '/')[[1]]
        allele_codes <- allele_codes[allele_codes != '.']
        allele_codes <- as.numeric(allele_codes)

        # Convert allele codes to sequence
        allele_code_key <- c(vcf_data$REF[index], strsplit(vcf_data$ALT[index], split = ',')[[1]])
        alleles <- allele_code_key[allele_codes + 1]

        # Return NA if not a simple SNP
        if (any(nchar(alleles) != 1)) {
            return(NA_character_)
        }

        # Convert alleles to a single sequence, using IUPAC ambiguity codes as needed
        return(iupac_key[paste0(sort(unique(alleles)), collapse = '')])
    })
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
