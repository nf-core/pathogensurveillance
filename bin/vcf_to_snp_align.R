#!/usr/bin/env Rscript

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- c(
    '/home/fosterz/projects/pathogensurveillance/work/cd/abe1a4b6a80d0f2d2598047d9c9104/_no_group_defined__GCF_000002765_6.vcffilter.vcf.gz',
    'deleteme.fasta'
)
names(args) <- c("vcf_path", "out_path")
args <- as.list(args)

# Read VCF file
header_data <- readLines(args$vcf_path, n = 1000)
header_line <- header_data[grepl(header_data, pattern = '^#CHROM')]
header <- strsplit(header_line, split = '\t')[[1]]
header[1] <- 'CHROM'
vcf_data <- read.delim(file = args$vcf_path, sep = '\t', comment.char = '#', header = FALSE)
colnames(vcf_data) <- header

# Make table with just alleles for each sample
sample_ids <- colnames(vcf_data)[10:length(colnames(vcf_data))]
gt_index <- vapply(strsplit(vcf_data$FORMAT, split = ':'), FUN.VALUE = integer(1), function(x) {
    which(x == 'GT')
})
vcf_data[sample_ids] <- lapply(sample_ids, function(samp_id) {
    vapply(seq_len(nrow(vcf_data)), FUN.VALUE = character(1), function(index) {
        strsplit(vcf_data[index, samp_id], split = ':')[[1]][gt_index[index]]
    })
})

# Convert allele numeric codes to AGCT + IUPAC ambiguity codes
iupac_key <- c(
    A="A",
    C="C",
    G="G",
    T="T",
    M="AC",
    R="AG",
    W="AT",
    S="CG",
    Y="CT",
    K="GT",
    V="ACG",
    H="ACT",
    D="AGT",
    B="CGT",
    N="ACGT"
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

# Remove variants with missing data
is_missing <- apply(vcf_data[sample_ids], MARGIN = 1, simplify = TRUE, function(x) any(is.na(x)))
vcf_data <- vcf_data[! is_missing, ]

# Remove indels relative to the reference
vcf_data <- vcf_data[nchar(vcf_data$REF) == 1, ]

# Convert to sequences
reference_seq <- paste0(vcf_data$REF, collapse = '')
names(reference_seq) <- "REF"
sample_seqs <- vapply(vcf_data[sample_ids], FUN.VALUE = character(1), paste0, collapse = '')
seqs <- c(reference_seq, sample_seqs)

# Write output FASTA file
writeLines(paste0('>', names(seqs), '\n', seqs), args$out_path)
