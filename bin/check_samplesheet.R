# This script takes 2 arguments:
#   1. The path to the CSV input to the pipeline supplied by the user
#   2. The path to the reformatted version of the CSV output by this script
#
# The first part of this script defines constants that might need to be changed in the future.

# Column names that can be used by the pipeline
# These will always be present and in this order in the output
# See the README.md for descriptions of each column
known_columns <- c(
    'sample_id',
    'sample_name',
    'shortread_1', 
    'shortread_2',
    'nanopore',
    'sra',
    'reference_id',
    'reference_name',
    'reference',
    'reference_refseq',
    'report_group',
    'color_by'
)

# Columns that must have a valid value in the input of this script
# For each vector in the list, at least one of the columns must have a value
required_input_columns <- list(
    c('shortread_1', 'shortread_2', 'nanopore', 'sra')
)

# Groups of columns in which only a single one should have a value. Regular expressions are allowed.
# If a regular expression matches multiple columns, then both matches can have a value
mutually_exclusive_columns <- list(
    c('reference', 'reference_refseq'),
    c('shortread_[12]', 'nanopore', 'sra')
)

# These are file extensions that are expected in the input data.
# These are (currently) only used to remove these extensions for file paths when they are used for IDs
known_extensions <- c(
    '.fastq',
    '.fastq.gz',
    '.fna',
    '.fna.gz',
    '.fasta',
    '.fasta.gz',
    '.fa',
    '.fa.gz'
)

# Regular expression for characters that cannot appear in IDs
invalid_id_char_pattern <- '[\\/:*?"<>| .]+'

# Name of default group if all samples do not have a group defined
defualt_group_full <- 'all'

# Name of default group if some samples do not have a group defined
defualt_group_partial <- '__other__'

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- as.list(args)
# args <- list('test/data/metadata_small.csv', 'test_out.csv') 
names(args) <- c('input_path', 'output_path')
metadata_original <- read.csv(args$input_path, check.names = FALSE)

# Check that there is data
if (nrow(metadata_original) == 0) {
    stop(call. = FALSE,
         'There are no rows in the input samplesheet.')
}

# Remove all whitespace
colnames(metadata_original) <- trimws(colnames(metadata_original))
metadata_original[] <- lapply(metadata_original, trimws)

# Check that required input columns are present
for (columns in required_input_columns) {
    if (! any(columns %in% colnames(metadata_original))) {
        stop(call. = FALSE,
             'At least one of the following columns must be present: ',
             paste0('"', columns, '"', collapse = ', '),
        )
    }
}

# Reorder columns and add any missing columns (initialized with NA)
empty_columns <- lapply(known_columns, function(column) {
    rep(NA_character_, nrow(metadata_original))
})
names(empty_columns) <- known_columns
metadata <- as.data.frame(empty_columns)
metadata[colnames(metadata_original)] <- metadata_original

# Check that required input columns have at least one value for each row
validate_required_input <- function(row_index, columns) {
    values <- unlist(metadata[row_index, columns])
    if (all(values == '' | is.na(values))) {
        stop(call. = FALSE, paste0(
            'At least one of the following columns must have a value on row ', row_index, ': ',
            paste0('"', columns, '"', collapse = ', ')
        ))
    }
}
for (row_index in 1:nrow(metadata)) {
    for (columns in required_input_columns) {
        validate_required_input(row_index, columns)
    }
}

# Validate mutually exclusive columns
validate_mutually_exclusive <- function(row_index, columns) {
    has_value <- unlist(lapply(columns, function(column) {
        col_pattern <- paste0('^', column, '$')
        matching_columns <- grep(colnames(metadata), pattern = col_pattern, value = TRUE)
        values <- metadata[row_index, matching_columns]
        return(any(values != '' & !is.na(values)))
    }))
    if (sum(has_value) > 1) {
        problem_columns <- columns[has_value]
        stop(call. = FALSE, paste0(
            'The following mutually exclusive columns all have values on row ', row_index, ': ',
            paste0('"', problem_columns, '"', collapse = ', '), '\n',
            'For this group of columns, only a single column type can have a value.'
        ))
    }
}
for (row_index in 1:nrow(metadata)) {
    for (columns in mutually_exclusive_columns) {
        validate_mutually_exclusive(row_index, columns)
    }
}

# Ensure sample IDs are present
shared_char <- function(col, end = FALSE) {
    col[col == ''] <- NA_character_
    reverse_char <- function(x) {
        split <- strsplit(x, split = "")
        reversed <- lapply(split, rev)
        return(unlist(lapply(reversed, paste, collapse = "")))
    }
    if (all(is.na(col))) {
        return("")
    }
    shorest_length <- min(unlist(lapply(col, nchar)), na.rm = TRUE)
    result <- ""
    if (end) {
        col <- reverse_char(col)
    }
    indexes <- 1:shorest_length
    for (index in indexes) {
        unique_starts <- unique(substr(col, start = 1, stop = index))
        if (length(unique_starts) == 1) {
            result <- unique_starts
        } else {
            break
        }
    }
    if (end) {
        result <- reverse_char(result)
    }
    return(result)
}
is_present <- function(x) {
    x != '' & ! is.na(x)
}
remove_shared <- function(x) {
    present_x <- x[is_present(x)]
    shared_start <- shared_char(present_x, end = FALSE)
    shared_end <- shared_char(present_x, end = TRUE)
    present_x <- sub(present_x, pattern = paste0('^', shared_start), replacement = '')
    present_x <- sub(present_x, pattern = paste0(shared_end, '$'), replacement = '')
    x[is_present(x)] <- present_x
    return(x)
}
remove_file_extensions <- function(x) {
    all_ext_pattern <- paste0('(', paste0(known_extensions, collapse = '|'), ')$')
    all_ext_pattern <- gsub(all_ext_pattern, pattern = '.', replacement = '\\.', fixed = TRUE)
    return(gsub(x, pattern = all_ext_pattern, replacement = ''))
}
shortread_ids <- unlist(lapply(1:nrow(metadata), function(row_index) {
    shortread_1 <- basename(metadata$shortread_1[row_index])
    shortread_2 <- basename(metadata$shortread_2[row_index])
    if (is_present(shortread_1) && is_present(shortread_2)) {
        remove_different_parts <- function(a, b) {
            a_split <- strsplit(shortread_1, split = '')[[1]]
            b_split <- strsplit(shortread_2, split = '')[[1]]
            paste0(a_split[a_split == b_split], collapse = '')
        }
        shortread <- remove_different_parts(shortread_1, shortread_2)
    } else if (is_present(shortread_1)) {
        shortread <- shortread_1
    } else if (is_present(shortread_2)) {
        shortread <- shortread_2
    } else {
        shortread <- NA_character_
    }
}))
id_sources <- list( # These are all possible sources of IDs, ordered by preference
    metadata$sample_id,
    metadata$sample_name,
    metadata$sra,
    remove_file_extensions(remove_shared(metadata$nanopore)),
    remove_file_extensions(remove_shared(shortread_ids))
)
metadata$sample_id <- unlist(lapply(1:nrow(metadata), function(row_index) { # Pick one replacement ID for each sample
    ids <- unlist(lapply(id_sources, `[`, row_index))
    return(ids[is_present(ids)][1])
}))

# Ensure sample names are present
metadata$sample_name <- unlist(lapply(1:nrow(metadata), function(row_index) {
    name <- metadata$sample_name[row_index]
    id <- metadata$sample_id[row_index]
    if (is_present(name)) {
        return(name)
    } else {
        return(id)
    }
}))

# Replace any characters in sample IDs that cannot be present in file names
metadata$sample_id <- gsub(metadata$sample_id, pattern = invalid_id_char_pattern, replacement = '_')

# Ensure sample IDs are unique
metadata$sample_id <- make.unique(metadata$sample_id, sep = '_')

# Ensure reference IDs are present
ref_id_sources <- list( # These are all possible sources of IDs, ordered by preference
    metadata$reference_id,
    metadata$reference_name,
    metadata$reference_refseq,
    remove_file_extensions(remove_shared(metadata$reference))
)
metadata$reference_id <- unlist(lapply(1:nrow(metadata), function(row_index) { # Pick one replacement ID for each sample
    ids <- unlist(lapply(ref_id_sources, `[`, row_index))
    return(ids[is_present(ids)][1])
}))

# Ensure reference names are present
metadata$reference_name <- unlist(lapply(1:nrow(metadata), function(row_index) {
    name <- metadata$reference_name[row_index]
    id <- metadata$reference_id[row_index]
    if (is_present(name)) {
        return(name)
    } else {
        return(id)
    }
}))

# Replace any characters in reference IDs that cannot be present in file names
metadata$reference_id <- gsub(metadata$reference_id, pattern = invalid_id_char_pattern, replacement = '_')

# Ensure reference IDs are unique
metadata$reference_id <- make.unique(metadata$reference_id, sep = '_')

# Add a default group for samples without a group defined
if (all(!is_present(metadata$report_group))) {
    metadata$report_group <- defualt_group_full
} else {
    metadata$report_group[!is_present(metadata$report_group)] <- defualt_group_partial
}

# Write output metadata
write.csv(metadata, file = args$output_path, row.names = FALSE, quote = FALSE, na = '')
