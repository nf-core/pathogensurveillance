#!/usr/bin/env Rscript

# This script takes 4 arguments:
#   1. The path to the per-sample CSV input to the pipeline supplied by the user
#   2. The path to the per-reference CSV input to the pipeline supplied by the user
#   3. The path to the reformatted version of the CSV per-sample output by this script
#   4. The path to the reformatted version of the CSV per-sample output by this script
#
# The first part of this script defines constants that might need to be changed in the future.

# Column names that can be used by the pipeline
# These will always be present and in this order in the output
# See the README.md for descriptions of each column
known_columns_samp <- c(
    'sample_id',
    'name',
    'description',
    'path',
    'path_2',
    'ncbi_accession',
    'ncbi_query',
    'ncbi_query_max',
    'sequence_type',
    'report_group_ids',
    'color_by',
    'ploidy',
    'enabled',
    'ref_group_ids',
    'ref_id',
    'ref_name',
    'ref_description',
    'ref_path',
    'ref_ncbi_accession',
    'ref_ncbi_query',
    'ref_ncbi_query_max',
    'ref_primary_usage',
    'ref_contextual_usage',
    'ref_color_by',
    'ref_enabled'
)
known_columns_ref <- c(
    'ref_group_ids',
    'ref_id',
    'ref_name',
    'ref_description',
    'ref_path',
    'ref_ncbi_accession',
    'ref_ncbi_query',
    'ref_ncbi_query_max',
    'ref_primary_usage',
    'ref_contextual_usage',
    'ref_color_by',
    'ref_enabled'
)

# Default values for columns
defaults_ref <- c(
    ref_ncbi_query_max = '30',
    ref_primary_usage = 'optional',
    ref_contextual_usage = 'optional',
    ref_enabled = TRUE
)
defaults_samp <- c(
    report_group_ids = '_no_group_defined_',
    enabled = TRUE,
    ploidy = '1',
    ncbi_query_max = '10',
    ref_ncbi_query_max = defaults_ref[['ref_ncbi_query_max']],
    ref_primary_usage = defaults_ref[['ref_primary_usage']],
    ref_contextual_usage = defaults_ref[['ref_contextual_usage']],
    ref_enabled = defaults_ref[['ref_enabled']]
)

# Columns that must have a valid value in the input of this script
# For each vector in the list, at least one of the columns must have a value
required_input_columns_samp <- list(
    c('path', 'path_2', 'ncbi_accession', 'ncbi_query')
)
required_input_columns_ref <- list(
    c('ref_path', 'ref_ncbi_accession', 'ref_ncbi_query')
)

# Groups of columns in which only a single one should have a value. Regular expressions are allowed.
# If a regular expression matches multiple columns, then both matches can have a value
mutually_exclusive_columns_samp <- list(
    c('ref_path', 'ref_ncbi_accession', 'ref_ncbi_query'),
    c('reads_?2?', 'ncbi_accession', 'ncbi_query')
)
mutually_exclusive_columns_ref <- list(
    c('ref_path', 'ref_ncbi_accession', 'ref_ncbi_query')
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
    '.fa.gz',
    '.fq',
    '.fq.gz'
)

# Types of sequencing supported by the pipeline.
# This is case insensitive
known_read_types <- c(
    'illumina',
    'nanopore',
    'pacbio'
)

# Regular expression for characters that cannot appear in IDs
invalid_id_char_pattern <- '[\\/:;*?"<>| .()-]+'

# Name of default group if all samples do not have a group defined
default_group_full <- 'all'

# Name of default report group if some samples do not have a group defined
default_group_partial <- '_no_group_defined_'

# Settings for how references can be used
valid_ref_usage_types <- c(
    'optional',
    'required',
    'excluded',
    'exclusive'
)

# Generally useful functions
is_present <- function(x) {
    x != '' & ! is.na(x)
}

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- as.list(args)
# args <- list('test/data/metadata/chaos_samples.csv', 'test/data/metadata/chaos_references.csv')
# args <- list('~/Downloads/sample_data_N273_14ncbigenomes.csv', '~/Downloads/ref_data.csv')
# args <- list('test/data/metadata/chaos_samples.csv')
metadata_original_samp <- read.csv(args[[1]], check.names = FALSE)
if (length(args) > 1) {
    metadata_original_ref <- read.csv(args[[2]], check.names = FALSE)
} else {
    metadata_original_ref <- data.frame(ref_path = character(0))
}
metadata_samp <- metadata_original_samp
metadata_ref <- metadata_original_ref

# Remove empty rows
remove_empty_rows <- function(metadata) {
    if (nrow(metadata) == 0) {
        return(metadata)
    }
    is_empty <- apply(metadata, MARGIN = 1, function(row) all(! is_present(row)))
    metadata[! is_empty, ]
}
metadata_samp <- remove_empty_rows(metadata_samp)
if (nrow(metadata_ref) > 0) {
    metadata_ref <- remove_empty_rows(metadata_ref)
}

# Check that there is data
if (nrow(metadata_samp) == 0) {
    stop(call. = FALSE,
         'There are no rows in the input samplesheet.')
}

# Validate column names
validate_col_names <- function(cols, known_columns) {
    # Trim whitespace
    modified_names <- trimws(cols)
    # Replace capital letters with lowercase 
    modified_names <- tolower(modified_names)
    # Replace spaces with underscores
    modified_names <- gsub(' +', '_', modified_names)
    # Add underscores in common missed locations
    underscore_replace_key <- known_columns
    names(underscore_replace_key) <- gsub ('_+', '', known_columns)
    colnames_to_replace <- underscore_replace_key[modified_names]
    modified_names[!is.na(colnames_to_replace)] <- colnames_to_replace[! is.na(colnames_to_replace)]
    # Apply changes to just known names
    cols[modified_names %in% known_columns] <- modified_names[modified_names %in% known_columns]
    return(cols)
}
colnames(metadata_samp) <- validate_col_names(colnames(metadata_samp), known_columns_samp)
if (nrow(metadata_ref) > 0) {
    colnames(metadata_ref) <- validate_col_names(colnames(metadata_ref), known_columns_ref)
}

# Remove empty columns
remove_empty_cols <- function(metadata, csv_name) {
    if (nrow(metadata) == 0) {
        return(metadata)
    }
    is_empty <- apply(metadata, MARGIN = 2, function(col) all(! is_present(col)))
    is_headerless <- ! is_present(colnames(metadata))
    if (any(is_headerless & ! is_empty)) {
        stop(call. = FALSE, paste0(
            'The following columns in the ', csv_name, ' input CSV have values but no header: ',
            paste0(which(is_headerless & ! is_empty), collapse = ', ')
        ))
    }
    metadata[, ! is_headerless]
}
metadata_samp <- remove_empty_cols(metadata_samp)
if (nrow(metadata_ref) > 0) {
    metadata_ref <- remove_empty_cols(metadata_ref)
}

# Remove all whitespace
remove_whitespace <- function(metadata) {
    metadata[] <- lapply(metadata, trimws)
    return(metadata)
}
metadata_samp <- remove_whitespace(metadata_samp)
if (nrow(metadata_ref) > 0) {
    metadata_ref <- remove_whitespace(metadata_ref)
}

# Replace NAs with empty stings
metadata_samp[] <- lapply(metadata_samp, function(x) {
    x[is.na(x)] <- ''
    return(x)
})
metadata_ref[] <- lapply(metadata_ref, function(x) {
    x[is.na(x)] <- ''
    return(x)
})

# Check that required input columns are present
check_required_cols <- function(metadata, required_cols, csv_name) {
    for (columns in required_cols) {
        if (! any(columns %in% colnames(metadata))) {
            stop(call. = FALSE,
                 'At least one of the following columns must be present in the ', csv_name, ' CSV: ',
                 paste0('"', columns, '"', collapse = ', ')
            )
        }
    }
}
check_required_cols(metadata_samp, required_input_columns_samp, 'sample data')
check_required_cols(metadata_ref, required_input_columns_ref, 'reference data')

# Check for duplicated columns
check_duplicated_cols <- function(metadata, known_cols, csv_name) {
    present_known_cols <- colnames(metadata)[colnames(metadata) %in% known_cols]
    duplicated_cols <- unique(present_known_cols[duplicated(present_known_cols)])
    if (length(duplicated_cols) > 0) {
        stop(call. = FALSE,
             'The following columns occur more than once in the ', csv_name, ' CSV: ',
             paste0('"', duplicated_cols, '"', collapse = ', ')
        )
    }
}
check_duplicated_cols(metadata_samp,  known_columns_samp, 'sample data')
check_duplicated_cols(metadata_ref,  known_columns_ref, 'reference data')

# Reorder columns and add any missing columns
reorder_and_add_cols <- function(metadata, known_columns) {
    empty_columns <- lapply(known_columns, function(column) {
        rep('', nrow(metadata))
    })
    names(empty_columns) <- known_columns
    output <- as.data.frame(empty_columns)
    output[colnames(metadata)] <- metadata
    return(output)
}
metadata_samp <- reorder_and_add_cols(metadata_samp, known_columns_samp)
metadata_ref <- reorder_and_add_cols(metadata_ref, known_columns_ref)

# Add default values for some columns
apply_defaults <- function(metadata, defaults) {
    metadata[names(defaults)] <- lapply(names(defaults), function(col) {
        ifelse(is_present(metadata[[col]]), metadata[[col]], defaults[col])
    })
    return(metadata)
}
metadata_samp <- apply_defaults(metadata_samp, defaults_samp)
if (nrow(metadata_ref) > 0) {
    metadata_ref <- apply_defaults(metadata_ref, defaults_ref)
}

# Validate mutually exclusive columns
validate_mutually_exclusive <- function(metadata, mutually_exclusive_columns, csv_name) {
    validate_mutually_exclusive_cols <- function(row_index, columns) {
        has_value <- unlist(lapply(columns, function(column) {
            col_pattern <- paste0('^', column, '$')
            matching_columns <- grep(colnames(metadata), pattern = col_pattern, value = TRUE)
            values <- metadata[row_index, matching_columns]
            return(any(values != '' & !is.na(values)))
        }))
        if (sum(has_value) > 1) {
            problem_columns <- columns[has_value]
            stop(call. = FALSE, paste0(
                'The following mutually exclusive columns in the ', csv_name, ' CSV all have values on row ', row_index, ': ',
                paste0('"', problem_columns, '"', collapse = ', '), '\n',
                'For this group of columns, only a single column type can have a value for each sample.'
            ))
        }
    }
    for (row_index in 1:nrow(metadata)) {
        for (columns in mutually_exclusive_columns_samp) {
            validate_mutually_exclusive_cols(row_index, columns)
        }
    }
}
validate_mutually_exclusive(metadata_samp, mutually_exclusive_columns_samp, 'sample data')
validate_mutually_exclusive(metadata_ref, mutually_exclusive_columns_ref, 'reference data')

# Remove any disabled rows
metadata_samp <- metadata_samp[as.logical(metadata_samp$enabled), ]
if (nrow(metadata_ref) > 0) {
    metadata_ref <- metadata_ref[as.logical(metadata_ref$ref_enabled), ]
}

# Validate color_by column and add back any original user-defined columns used
validate_color_by <- function(metadata, color_by_col, known_cols, csv_name, sep = ';') {
    split_color_by <- strsplit(metadata[[color_by_col]], split = sep)
    split_color_by <- lapply(split_color_by, trimws)
    all_color_by_cols <- unique(unlist(split_color_by))
    missing_cols <- all_color_by_cols[! all_color_by_cols %in% colnames(metadata)]
    if (length(missing_cols) > 0) {
        stop(call. = FALSE, paste0(
            'The following columns in the ', csv_name, ' CSV referenced by the "', color_by_col, 
            '" column do not exist: ', paste0('"', missing_cols, '"', collapse = ', '), '\n'
        ))
    }
    return(unlist(lapply(split_color_by, paste0, collapse = sep)))
}
metadata_samp$color_by <- validate_color_by(metadata_samp, 'color_by', known_columns_samp, 'sample data')
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_color_by <- validate_color_by(metadata_ref, 'ref_color_by', known_columns_ref, 'reference data')
}

# Move reference data from the sample metadata to the reference metadata
ref_in_samp_data <- metadata_samp[, known_columns_ref]
has_ref_data <- apply(ref_in_samp_data[, unlist(required_input_columns_ref)], MARGIN = 1, function(x) any(is_present(x)))
ref_data_addition <- ref_in_samp_data[has_ref_data, ]
ref_data_addition <- ref_data_addition[, colnames(metadata_ref)]
metadata_ref <- rbind(metadata_ref, ref_data_addition)

# Validate usage columns 
validate_usage_col <- function(metadata, col) {
    unlist(lapply(1:nrow(metadata), function(index) {
        value <- tolower(trimws(metadata[[col]][index]))
        if (! value %in% valid_ref_usage_types) {
            stop(call. = FALSE, paste0(
                'The value "', value, '" on row ', index, ' column "', col, '" is not valid. It must be one of "',
                paste0(valid_ref_usage_types, collapse = '", "'), '"'
            ))
        }
        return(value)
    }))
 }
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_primary_usage <- validate_usage_col(metadata_ref, 'ref_primary_usage')
    metadata_ref$ref_contextual_usage <- validate_usage_col(metadata_ref, 'ref_contextual_usage')
}

# Validate ploidy column
for (index in 1:nrow(metadata_samp)) {
    if (! grepl(metadata_samp$ploidy[index], pattern = '^[0-9]+$')) {
        stop(call. = FALSE, paste0(
            'The value for the ploidy column "', metadata_samp$ploidy[index], 
            '" of the sample metadata CSV on row ', index, ' is not numeric.'
        ))
    }
}

# Convert NCBI sample queries to a list of SRA run accessions
get_ncbi_sra_runs <- function(query) {
    if (query == '') {
        return(NULL)
    }
    search_result <- rentrez::entrez_search(db = 'sra', query, retmax = 300)
    summary_result <- rentrez::entrez_summary(db = 'sra', search_result$ids, retmax = 300)
    if (length(search_result$ids) == 1) {
        summary_result <- list(summary_result)
    }
    run_ids <- unlist(lapply(summary_result, function(x) {
        if (length(x$runs) > 1) {
            warning('The SRA accession ', x$uid, ' has multiple runs associated with it. Only the first will be used.')
        }
        run_xml <- x$runs[1]
        gsub(run_xml, pattern = '.+ acc="(.+?)" .+', replacement = '\\1')
    }))
    sequence_instrument <- unlist(lapply(summary_result, function(x) {
        gsub(x$expxml[1], pattern = '.+<Platform instrument_model="(.+?)">(.+?)</Platform>.+', replacement = '\\1')
    }))
    sequence_type <- unlist(lapply(summary_result, function(x) {
        gsub(x$expxml[1], pattern = '.+<Platform instrument_model="(.+?)">(.+?)</Platform>.+', replacement = '\\2')
    }))
    name <- unlist(lapply(summary_result, function(x) {
        gsub(x$expxml[1], pattern = '.+ ScientificName="(.+?)"/>.+', replacement = '\\1')
    }))
    description <- unlist(lapply(summary_result, function(x) {
        gsub(x$expxml[1], pattern = '.+<Title>(.+?)</Title>.+', replacement = '\\1')
    }))
    output <- data.frame(
        ncbi_accession = run_ids,
        sequence_type = sequence_type,
        name = name,
        description = paste0(description, ' (', run_ids, ')')
    )
    rownames(output) <- NULL
    return(output)
}
unique_queries <- unique(metadata_samp$ncbi_query)
unique_queries <- unique_queries[unique_queries != '']
ncbi_result <- lapply(unique_queries, get_ncbi_sra_runs)
names(ncbi_result) <- unique_queries
new_sample_data <- do.call(rbind, lapply(which(is_present(metadata_samp$ncbi_query)), function(index) {
    query_data <- ncbi_result[[metadata_samp$ncbi_query[index]]]
    query_max <- metadata_samp$ncbi_query_max[index]
    if (endsWith(query_max, '%')) {
        query_max_prop <- as.numeric(gsub(query_max, pattern = '%$', replacement = '')) / 100
        query_max <- max(c(1, nrow(query_data) * query_max_prop))
    } else {
        query_max <- min(c(nrow(query_data), as.numeric(query_max)))
    }
    query_data <- query_data[sample(1:nrow(query_data), query_max), ]
    output <- metadata_samp[rep(index, nrow(query_data)), ]
    output$sample_id <- paste0(metadata_samp$sample_id[index], query_data$ncbi_accession)
    output$name <- paste0(metadata_samp$name[index], query_data$name)
    output$description <- paste0(metadata_samp$description[index], query_data$description)
    output$ncbi_accession <- query_data$ncbi_accession
    output$sequence_type <- query_data$sequence_type
    rownames(output) <- NULL
    output
}))
metadata_samp <- rbind(
    metadata_samp[! is_present(metadata_samp$ncbi_query), ], 
    new_sample_data
)

# Convert NCBI reference queries to a list of assembly accessions
get_ncbi_genomes <- function(query) {
    search_result <- rentrez::entrez_search(db = 'assembly', query, retmax = 300)
    summary_result <- rentrez::entrez_summary(db = 'assembly', search_result$ids, retmax = 300)
    if (length(search_result$ids) == 1) {
        summary_result <- list(summary_result)
    }
    acc_ids <- unlist(lapply(summary_result, function(x) {
        run_xml <- x$assemblyaccession[1]
        gsub(run_xml, pattern = '.+ acc="(.+?)" .+', replacement = '\\1')
    }))
    output <- data.frame(
        ref_id = unlist(lapply(summary_result, function(x) x$assemblyaccession)),
        ref_name = unlist(lapply(summary_result, function(x) paste0(x$speciesname, ' ', x$assemblyname))),
        ref_description = unlist(lapply(summary_result, function(x) paste0(x$organism, ' ', x$assemblyname, ' (', x$assemblyaccession, ')'))),
        ref_ncbi_accession = unlist(lapply(summary_result, function(x) x$assemblyaccession))
    )
    rownames(output) <- NULL
    return(output)
}
unique_queries <- unique(metadata_ref$ref_ncbi_query)
unique_queries <- unique_queries[unique_queries != '']
ncbi_result <- lapply(unique_queries, get_ncbi_genomes)
names(ncbi_result) <- unique_queries
new_ref_data <- do.call(rbind, lapply(which(is_present(metadata_ref$ref_ncbi_query)), function(index) {
    query_data <- ncbi_result[[metadata_ref$ref_ncbi_query[index]]]
    query_max <- metadata_ref$ref_ncbi_query_max[index]
    if (endsWith(query_max, '%')) {
        query_max_prop <- as.numeric(gsub(query_max, pattern = '%$', replacement = '')) / 100
        query_max <- max(c(1, nrow(query_data) * query_max_prop))
    } else {
        query_max <- min(c(nrow(query_data), as.numeric(query_max)))
    }
    query_data <- query_data[sample(1:nrow(query_data), query_max), ]
    output <- metadata_ref[rep(index, nrow(query_data)), ]
    output$ref_id <- paste0(metadata_ref$ref_id[index], query_data$ref_ncbi_accession)
    output$ref_name <- paste0(metadata_ref$ref_name[index], query_data$ref_name)
    output$ref_description <- paste0(metadata_ref$ref_description[index], query_data$ref_description)
    output$ref_ncbi_accession <- query_data$ref_ncbi_accession
    rownames(output) <- NULL
    return(output)
}))
metadata_ref <- rbind(
    metadata_ref[! is_present(metadata_ref$ref_ncbi_query), ], 
    new_ref_data
)

# Check that required input columns have at least one value for each row
validate_required_input <- function(metadata, required_input_columns, csv_name) {
    validate_required_input_cols <- function(row_index, columns) {
        values <- unlist(metadata[row_index, columns])
        if (all(values == '' | is.na(values))) {
            stop(call. = FALSE, paste0(
                'At least one of the following columns in the ', csv_name, ' CSV must have a value on row ', row_index, ': ',
                paste0('"', columns, '"', collapse = ', ')
            ))
        }
    }
    for (row_index in seq_along(rownames(metadata))) {
        for (columns in required_input_columns) {
            validate_required_input_cols(row_index, columns)
        }
    }
}
validate_required_input(metadata_samp, required_input_columns_samp, 'sample data')
if (nrow(metadata_ref) > 0) {
    validate_required_input(metadata_ref, required_input_columns_ref, 'reference data')
}

# Ensure sample/reference IDs are present
shared_char <- function(col, end = FALSE) {
    col[col == ''] <- NA_character_
    reverse_char <- function(x) {
        split <- strsplit(x, split = "")
        reversed <- lapply(split, rev)
        return(unlist(lapply(reversed, paste, collapse = "")))
    }
    if (all(is.na(col))) {
        return('')
    }
    shorest_length <- min(unlist(lapply(col, nchar)), na.rm = TRUE)
    result <- ''
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

remove_shared <- function(x) {
    if (length(x) > 1) {
        present_x <- x[is_present(x)]
        shared_start <- shared_char(present_x, end = FALSE)
        shared_end <- shared_char(present_x, end = TRUE)
        present_x <- sub(present_x, pattern = paste0('^', shared_start), replacement = '')
        present_x <- sub(present_x, pattern = paste0(shared_end, '$'), replacement = '')
        x[is_present(x)] <- present_x
    } else {
        x = basename(x)
    }
    return(x)
}
remove_file_extensions <- function(x) {
    all_ext_pattern <- paste0('(', paste0(known_extensions, collapse = '|'), ')$')
    all_ext_pattern <- gsub(all_ext_pattern, pattern = '.', replacement = '\\.', fixed = TRUE)
    return(gsub(x, pattern = all_ext_pattern, replacement = ''))
}

reads_ids <- unlist(lapply(1:nrow(metadata_samp), function(row_index) {
    reads_1 <- basename(metadata_samp$path[row_index])
    reads_2 <- basename(metadata_samp$path_2[row_index])
    if (is_present(reads_1) && is_present(reads_2)) {
        remove_different_parts <- function(a, b) {
            a_split <- strsplit(reads_1, split = '')[[1]]
            b_split <- strsplit(reads_2, split = '')[[1]]
            paste0(a_split[a_split == b_split], collapse = '')
        }
        shortread <- remove_different_parts(reads_1, reads_2)
    } else if (is_present(reads_1)) {
        shortread <- reads_1
    } else if (is_present(reads_2)) {
        shortread <- reads_2
    } else {
        shortread <- ''
    }
}))

id_sources_samp <- list( # These are all possible sources of IDs, ordered by preference
    metadata_samp$sample_id,
    metadata_samp$sample_name,
    metadata_samp$ncbi_accession,
    remove_file_extensions(remove_shared(reads_ids))
)
metadata_samp$sample_id <- unlist(lapply(1:nrow(metadata_samp), function(row_index) { # Pick one replacement ID for each sample
    ids <- unlist(lapply(id_sources_samp, `[`, row_index))
    return(ids[is_present(ids)][1])
}))

id_sources_ref <- list( # These are all possible sources of IDs, ordered by preference
    metadata_ref$ref_id,
    metadata_ref$ref_name,
    metadata_ref$ref_ncbi_accession,
    metadata_ref$ref_path
)
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_id <- unlist(lapply(1:nrow(metadata_ref), function(row_index) { # Pick one replacement ID for each sample
        ids <- unlist(lapply(id_sources_ref, `[`, row_index))
        return(ids[is_present(ids)][1])
    }))
}

# Ensure sample/reference names and descriptions are present
ensure_sample_names <- function(names, ids) {
    unlist(lapply(seq_along(names), function(index) {
        if (is_present(names[index])) {
            return(names[index])
        } else {
            return(ids[index])
        }
    }))
}
metadata_samp$name <- ensure_sample_names(metadata_samp$name, metadata_samp$sample_id)
metadata_samp$description <- ensure_sample_names(metadata_samp$description, metadata_samp$name)
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_name <- ensure_sample_names(metadata_ref$ref_name, metadata_ref$ref_id)
}
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_description <- ensure_sample_names(metadata_ref$ref_description, metadata_ref$ref_name)
}

# Replace any characters in IDs that cannot be present in file names
make_ids_ok_for_file_names <- function(ids) {
    trimws(gsub(ids, pattern = invalid_id_char_pattern, replacement = '_'))
}
metadata_samp$sample_id <- make_ids_ok_for_file_names(metadata_samp$sample_id)
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_id <- make_ids_ok_for_file_names(metadata_ref$ref_id)
}

# Ensure that the same sample ID is not used for different sets of data
make_ids_unique <- function(metadata, id_col, other_cols) {
    # Find which IDs need to be changed
    subset <- metadata[, c(id_col, other_cols)]
    subset$row_num <- seq_along(rownames(subset))
    unique_ids <- unique(subset[[id_col]])
    unique_ids <- unique_ids[is_present(unique_ids)]
    id_key <- lapply(unique_ids, function(id) {
        same_id_rows <- subset[subset[[id_col]] == id, ]
        split_by_other <- split(same_id_rows, apply(same_id_rows[, other_cols], 1, paste0, collapse = ''))
        names(split_by_other) <- make.unique(rep(id, length(split_by_other)), sep = '_')
        new_id_key <- lapply(names(split_by_other), function(new_id) { # Get list of data.frames with new IDs and associated row numbers
            data.frame(new_id = new_id, row_num = split_by_other[[new_id]]$row_num)
        })
        new_id_key <- do.call(rbind, new_id_key) # combine list of data.frames to a single one
        return(new_id_key)
    })
    id_key <- do.call(rbind, id_key) # combine list of data.frames to a single one
    # Apply the changes
    metadata[id_key$row_num, id_col] <- id_key$new_id
    return(metadata)
}
metadata_samp <- make_ids_unique(metadata_samp, id_col = 'sample_id', other_cols = c('path', 'path_2', 'ncbi_accession'))
if (nrow(metadata_ref) > 0) {
    metadata_ref <- make_ids_unique(metadata_ref, id_col = 'ref_id', other_cols = c('ref_path', 'ref_ncbi_accession'))
}

# Ensure references and samples do not share ids
is_shared <- metadata_ref$ref_id %in% metadata_samp$sample_id
metadata_ref$ref_id[is_shared] <- paste0(metadata_ref$ref_id[is_shared], '_ref')

# Add a default report group for samples without a group defined
if (all(!is_present(metadata_samp$report_group_ids))) {
    metadata_samp$report_group_ids <- default_group_full
} else {
    metadata_samp$report_group_ids[!is_present(metadata_samp$report_group_ids)] <- default_group_partial
}

# Add reference id to the reference group ids, so that references can be referred to by id even if a group id is not defined.
metadata_ref$ref_group_ids <- ifelse(
    metadata_ref$ref_group_ids == '',
    metadata_ref$ref_id,
    paste(metadata_ref$ref_group_ids, metadata_ref$ref_id, sep = ';')
)

# Make report/reference group ids usable as file names
make_group_ids_ok_for_file_names <- function(group_ids, sep = ';') {
    unlist(lapply(group_ids, function(ids) {
        split_ids <- strsplit(ids, split = sep)[[1]]
        split_ids <- make_ids_ok_for_file_names(split_ids)
        paste(split_ids, collapse = sep)
    }))
}
metadata_samp$report_group_ids <- make_group_ids_ok_for_file_names(metadata_samp$report_group_ids)
metadata_samp$ref_group_ids <- make_group_ids_ok_for_file_names(metadata_samp$ref_group_ids)
if (nrow(metadata_ref) > 0) {
    metadata_ref$ref_group_ids <- make_group_ids_ok_for_file_names(metadata_ref$ref_group_ids)
}

# Check that reference groups in sample metadata are present in the reference metadata
if (nrow(metadata_ref)) {
    all_ref_group_ids <- unique(unlist(strsplit(metadata_ref$ref_group_ids, split = ';')))
    for (index in 1:nrow(metadata_samp)) {
        split_ids <- strsplit(metadata_samp$ref_group_ids[index], split = ';')[[1]]
        invalid_ids <- split_ids[! split_ids %in% all_ref_group_ids]
        if (length(invalid_ids) > 0) {
            stop(call. = FALSE, paste0(
                'The reference group ID "', invalid_ids[1], '" used in row ', index, ' in the sample metadata CSV',
                ' is not defined in the reference metadata CSV. All values in the "ref_group_ids" column in the',
                ' sample metadata CSV must be present in the "ref_group_ids" or "ref_id" columns of the reference',
                ' metadata CSV.'
            ))
        }
    }
}

# Look up sequence type for NCBI accessions without a sequence type defined
lookup_sequence_type <- function(id) {
    search_result <- rentrez::entrez_search(db = 'sra', id)
    summary_result <- rentrez::entrez_summary(db = 'sra', search_result$ids)
    if (length(search_result$ids) != 1) {
        stop(call. = FALSE, paste0(
            'Could not look up sequence type for id "', id, '". Specify it in the sample metadata CSV.'
        ))
    }
    gsub(summary_result$expxml, pattern = '.+<Platform instrument_model="(.+?)">(.+?)</Platform>.+', replacement = '\\2')
}
undefined_accessions <- unique(metadata_samp$ncbi_accession[! is_present(metadata_samp$sequence_type) & is_present(metadata_samp$ncbi_accession)])
type_replace_key <- lapply(undefined_accessions, lookup_sequence_type)
names(type_replace_key) <- undefined_accessions
is_undefined <- metadata_samp$ncbi_accession %in% undefined_accessions
metadata_samp$sequence_type[is_undefined] <- type_replace_key[metadata_samp$ncbi_accession[is_undefined]]

# Validate sequence_type column
metadata_samp$sequence_type <- unlist(lapply(seq_along(metadata_samp$sequence_type), function(index) {
    is_seq_type <- unlist(lapply(known_read_types, function(type) {
        grepl(metadata_samp$sequence_type[index], pattern = type, ignore.case = TRUE)
    }))
    if (sum(is_seq_type) == 0) {
        stop(call. = FALSE, paste0(
            'The value in the "sequence_type" column on row ', index, ' does not contain a known sequence type. ',
            'One of the following words must appear (case insensitive):\n',
            paste0('"', known_read_types, '"', collapse = ', '), '\n'
        ))
    }
    if (sum(is_seq_type) > 1) {
        stop(call. = FALSE, paste0(
            'The value in the "sequence_type" column on row ', index, ' contains the names of multiple sequencing types. ',
            'Exactly one of the following words must appear (case insensitive):\n',
            paste0('"', known_read_types, '"', collapse = ', '), '\n'
        ))
    }
    return(known_read_types[is_seq_type])
}))

# Add row for each group for samples/references with multiple groups
duplicate_rows_by_id_list <- function(metadata, id_col) {
    group_ids <- strsplit(metadata[[id_col]], split = ';')
    n_group_ids <- unlist(lapply(group_ids, length))
    group_ids[n_group_ids == 0] <- ''
    n_group_ids <- unlist(lapply(group_ids, length))
    metadata <- metadata[rep(1:nrow(metadata), n_group_ids), ]
    metadata[[id_col]] <- unlist(group_ids)
    rownames(metadata) <- NULL
    return(metadata)
}
metadata_samp <- duplicate_rows_by_id_list(metadata_samp, 'report_group_ids')
if (nrow(metadata_ref) > 0) {
    metadata_ref <- duplicate_rows_by_id_list(metadata_ref, 'ref_group_ids')
}

# Convert reference groups to reference ids in the sample data
metadata_samp$ref_ids <- unlist(lapply(metadata_samp$ref_group_ids, function(group_ids) {
    ref_ids <- metadata_ref$ref_id[metadata_ref$ref_group_ids %in% strsplit(group_ids, split = ';')]
    paste(ref_ids, collapse = ';')
}))

# Remove unneeded columns
metadata_samp$enabled <- NULL
metadata_samp$ref_group_ids <- NULL
metadata_samp$ref_id <- NULL
metadata_samp$ref_name <- NULL
metadata_samp$ref_description <- NULL
metadata_samp$ref_path <- NULL
metadata_samp$ref_ncbi_accession <- NULL
metadata_samp$ref_ncbi_query <- NULL
metadata_samp$ref_ncbi_query_max <- NULL
metadata_samp$ref_primary_usage <- NULL
metadata_samp$ref_color_by <- NULL
metadata_samp$ref_enabled <- NULL
metadata_ref$ref_group_ids <- NULL
metadata_ref$ref_enabled <- NULL

# Make rows unique
metadata_samp <- unique(metadata_samp)
if (nrow(metadata_ref) > 0) {
    metadata_ref <- unique(metadata_ref)
}

# Write output metadata
write.csv(metadata_samp, file = 'sample_metadata.csv', row.names = FALSE, na = '')
write.csv(metadata_ref, file = 'reference_metadata.csv', row.names = FALSE, na = '')
