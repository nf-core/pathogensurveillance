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
    'sequence_type',
    'report_group_ids',
    'color_by',
    'enabled',
    'ploidy',
    'ref_id',
    'ref_name',
    'ref_description',
    'ref_group_ids',
    'ref_path',
    'ref_ncbi_accession',
    'ref_ncbi_query',
    'ref_primary_usage',
    'ref_contextual_usage',
    'ref_color_by',
    'ref_enabled'
)
known_columns_ref <- c(
    'ref_id',
    'ref_name',
    'ref_description',
    'ref_group_ids',
    'ref_path',
    'ref_ncbi_accession',
    'ref_ncbi_query',
    'ref_primary_usage',
    'ref_contextual_usage',
    'ref_color_by',
    'ref_enabled'
)

# Columns that must have a valid value in the input of this script
# For each vector in the list, at least one of the columns must have a value
required_input_columns_samp <- list(
    c('path', 'path_2', 'ncbi_accession', 'ncbi_query'),
    c('sequence_type')
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
invalid_id_char_pattern <- '[\\/:*?"<>| .-]+'

# Name of default group if all samples do not have a group defined
default_group_full <- 'all'

# Name of default group if some samples do not have a group defined
default_group_partial <- '__other__'

# Prefix added to column names to distinguish modified columns from user-supplied columns
user_column_name_prefix <- '__user_'

# Parse inputs
args <- commandArgs(trailingOnly = TRUE)
args <- as.list(args)
args <- list('test/data/metadata/chaos_samples.csv', 'test/data/metadata/chaos_references.csv', 'test_out_samp.csv',  'test_out_ref.csv')
names(args) <- c('input_path_samp', 'input_path_ref', 'output_path_samp', 'output_path_ref')
metadata_original_samp <- read.csv(args$input_path_samp, check.names = FALSE)
metadata_original_ref <- read.csv(args$input_path_ref, check.names = FALSE)

# Remove empty rows
remove_empty_rows <- function(metadata) {
    is_empty <- apply(metadata, MARGIN = 1, function(row) all(is.na(row) | row == ''))
    metadata[! is_empty, ]
}
metadata_original_samp <- remove_empty_rows(metadata_original_samp)
metadata_original_ref <- remove_empty_rows(metadata_original_ref)

# Remove empty columns
remove_empty_cols <- function(metadata) {
    is_empty <- apply(metadata, MARGIN = 2, function(row) all(is.na(row) | row == ''))
    metadata[, ! is_empty]
}
metadata_original_samp <- remove_empty_cols(metadata_original_samp)
metadata_original_ref <- remove_empty_cols(metadata_original_ref)

# Check that there is data
if (nrow(metadata_original_samp) == 0) {
    stop(call. = FALSE,
         'There are no rows in the input samplesheet.')
}

# Remove all whitespace
remove_whitespace <- function(metadata) {
    colnames(metadata) <- trimws(colnames(metadata))
    metadata[] <- lapply(metadata, trimws)
    return(metadata)
}
metadata_original_samp <- remove_whitespace(metadata_original_samp)
metadata_original_ref <- remove_whitespace(metadata_original_ref)

# Preserve original column names
unmodified_data_samp <- metadata_original_samp
unmodified_data_ref <- metadata_original_ref

# Replace capital letters with lowercase in colnames
colnames(metadata_original_samp) <- tolower(colnames(metadata_original_samp))
colnames(metadata_original_ref) <- tolower(colnames(metadata_original_ref))

# Replace spaces with underscores in colnames
colnames(metadata_original_samp) <- gsub(' +', '_', colnames(metadata_original_samp))
colnames(metadata_original_ref) <- gsub(' +', '_', colnames(metadata_original_ref))

# Add underscores in common missed locations
add_underscores <- function(metadata, known_columns) {
    underscore_replace_key <- known_columns
    names(underscore_replace_key) <- gsub ('_+', '', known_columns)
    colnames_to_replace <- underscore_replace_key[colnames(metadata)]
    colnames(metadata)[!is.na(colnames_to_replace)] <- colnames_to_replace[! is.na(colnames_to_replace)]
    return(metadata)
}
metadata_original_samp <- add_underscores(metadata_original_samp, known_columns_samp)
metadata_original_ref <- add_underscores(metadata_original_ref, known_columns_ref)

# Replace NAs with empty stings
metadata_original_samp[] <- lapply(metadata_original_samp, function(x) {
    x[is.na(x)] <- ''
    return(x)
})
metadata_original_ref[] <- lapply(metadata_original_ref, function(x) {
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
check_required_cols(metadata_original_samp, required_input_columns_samp, 'sample data')
check_required_cols(metadata_original_ref, required_input_columns_ref, 'reference data')

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
check_duplicated_cols(metadata_original_samp,  known_columns_samp, 'sample data')
check_duplicated_cols(metadata_original_ref,  known_columns_ref, 'reference data')

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
metadata_samp <- reorder_and_add_cols(metadata_original_samp, known_columns_samp)
metadata_ref <- reorder_and_add_cols(metadata_original_ref, known_columns_ref)

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

# Convert NCBI sample queries to a list of SRA run accessions
get_ncbi_sra_runs <- function(query) {
    if (query == '') {
        return(NULL)
    }
    search_result <- rentrez::entrez_search(db = 'sra', query)
    summary_result <- rentrez::entrez_summary(db = 'sra', search_result$ids)
    if (length(search_result$ids) == 1) {
        summary_result <- list(summary_result)
    }
    run_ids <- unlist(lapply(summary_result, function(x) {
        if (length(x$runs) > 1) {
            warning('The SRA accession ', x$uid, ' has multiple runs associated with it. Only the first will be used.')
        }
        run_xml <- x$runs[1]
        gsub(run_xml, pattern = '.+ acc="([a-zA-Z0-9.]+)" .+', replacement = '\\1')
    }))
    sequence_instrument <- unlist(lapply(summary_result, function(x) {
        gsub(x$expxml[1], pattern = '.+<Platform instrument_model="([a-zA-Z0-9. ]+)">([a-zA-Z0-9.]+)</Platform>.+', replacement = '\\1')
    }))
    sequence_type <- unlist(lapply(summary_result, function(x) {
        gsub(x$expxml[1], pattern = '.+<Platform instrument_model="([a-zA-Z0-9. ]+)">([a-zA-Z0-9.]+)</Platform>.+', replacement = '\\2')
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

new_sample_data <- do.call(rbind, lapply(which(metadata_samp$ncbi_query != ""), function(index) {
    query_data <- ncbi_result[[metadata_samp$ncbi_query[index]]]
    query_data$sample_id <- paste0(metadata_samp$sample_id[index], query_data$ncbi_accession)
    query_data$name <- paste0(metadata_samp$name[index], query_data$name)
    query_data$description <- paste0(metadata_samp$description[index], query_data$description)
    query_data
}))

format_new_data_as_old <- function(new_metadata, old_metadata) {
    empty_columns <- lapply(colnames(old_metadata), function(column) {
        rep('', nrow(new_metadata))
    })
    names(empty_columns) <- colnames(old_metadata)
    output <- as.data.frame(empty_columns)
    output[colnames(new_metadata)] <- new_metadata
    return(output)
}

metadata_samp <- rbind(
    metadata_samp[metadata_samp$ncbi_query == '', ], 
    format_new_data_as_old(new_sample_data, metadata_samp)
)



# Convert NCBI reference queries to a list of assembly accessions
get_ncbi_genomes <- function(query) {
    search_result <- rentrez::entrez_search(db = 'assembly', query)
    summary_result <- rentrez::entrez_summary(db = 'assembly', search_result$ids)
    if (length(search_result$ids) == 1) {
        summary_result <- list(summary_result)
    }
    acc_ids <- unlist(lapply(summary_result, function(x) {
        run_xml <- x$assemblyaccession[1]
        gsub(run_xml, pattern = '.+ acc="([a-zA-Z0-9.]+)" .+', replacement = '\\1')
    }))
    output <- data.frame(
        ref_id = unlist(lapply(summary_result, function(x) x$assemblyaccession)),
        ref_name = unlist(lapply(summary_result, function(x) paste0(x$speciesname, x$assemblyname))),
        ref_description = unlist(lapply(summary_result, function(x) paste0(x$organism, ' ', x$assemblyname, ' (', x$assemblyaccession, ')'))),
        ncbi_accession = unlist(lapply(summary_result, function(x) x$assemblyaccession))
    )
    rownames(output) <- NULL
    return(output)
}

unique_queries <- unique(metadata_ref$ref_ncbi_query)
unique_queries <- unique_queries[unique_queries != '']
ncbi_result <- lapply(unique_queries, get_ncbi_genomes)
names(ncbi_result) <- unique_queries

new_ref_data <- do.call(rbind, lapply(which(metadata_ref$ref_ncbi_query != ""), function(index) {
    query_data <- ncbi_result[[metadata_ref$ref_ncbi_query[index]]]
    query_data$sample_id <- paste0(metadata_ref$sample_id[index], query_data$ncbi_accession)
    query_data$name <- paste0(metadata_ref$ref_name[index], query_data$ref_name)
    query_data$description <- paste0(metadata_ref$ref_description[index], query_data$ref_description)
    query_data
}))


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
    for (row_index in 1:nrow(metadata)) {
        for (columns in required_input_columns) {
            validate_required_input_cols(row_index, columns)
        }
    }
}
validate_required_input(metadata_original_samp, required_input_columns_samp, 'sample data')
validate_required_input(metadata_original_ref, required_input_columns_ref, 'reference data')


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
is_present <- function(x) {
    x != '' & ! is.na(x)
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

reads_ids <- unlist(lapply(1:nrow(metadata), function(row_index) {
    reads_1 <- basename(metadata$reads_1[row_index])
    reads_2 <- basename(metadata$reads_2[row_index])
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
id_sources <- list( # These are all possible sources of IDs, ordered by preference
    metadata$sample_id,
    metadata$sample_name,
    metadata$reads_sra,
    remove_file_extensions(remove_shared(reads_ids))
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

# Ensure that the same sample ID is not used for different sets of data
make_ids_unique <- function(metadata, id_col, other_cols) {
    # Find which IDs need to be changed
    subset <- metadata[, c(id_col, other_cols)]
    subset$row_num <- 1:nrow(subset)
    unique_ids <- unique(subset[[id_col]])
    unique_ids <- unique_ids[is_present(unique_ids)]
    id_key <- lapply(unique_ids, function(id) {
        print(id)
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
metadata <- make_ids_unique(metadata, id_col = 'sample_id', other_cols = c('reads_1', 'reads_2', 'reads_sra'))

# Ensure reference IDs are present
ref_id_sources <- list( # These are all possible sources of IDs, ordered by preference
    metadata$reference_id,
    metadata$reference_name,
    metadata$reference_refseq,
    remove_file_extensions(remove_shared(metadata$reference))
)
metadata$reference_id <- unlist(lapply(1:nrow(metadata), function(row_index) { # Pick one replacement ID for each sample
    ids <- unlist(lapply(ref_id_sources, `[`, row_index))
    id <- ids[is_present(ids)][1]
    if (is.na(id)) {
        id <- ''
    }
    return(id)
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

# Remove reference names and IDs for samples that have no reference specified
no_reference_data <- ! is_present(metadata$reference) & ! is_present(metadata$reference_refseq)
metadata$reference_id[no_reference_data] <- ''
metadata$reference_name[no_reference_data] <- ''

# Replace any characters in reference IDs that cannot be present in file names
metadata$reference_id <- gsub(metadata$reference_id, pattern = invalid_id_char_pattern, replacement = '_')

# Ensure that the same reference ID is not used for different sets of data
metadata <- make_ids_unique(metadata, id_col = 'reference_id', other_cols = c('reference', 'reference_refseq'))

# Ensure references and samples do not share ids
is_shared <- metadata$reference_id %in% metadata$sample_id
metadata$reference_id[is_shared] <- paste0(metadata$reference_id[is_shared], '_ref')

# Add a default group for samples without a group defined
if (all(!is_present(metadata$report_group))) {
    metadata$report_group <- default_group_full
} else {
    metadata$report_group[!is_present(metadata$report_group)] <- default_group_partial
}

# Remove whitespace in report group ids
metadata$report_group <- trimws(metadata$report_group)
metadata$report_group <- gsub(metadata$report_group, pattern = '[[:space:]]+;[[:space:]]+', replacement = ';')

# Validate reads_type column
metadata$reads_type <- unlist(lapply(seq_along(metadata$reads_type), function(index) {
    is_seq_type <- unlist(lapply(known_read_types, function(type) {
        grepl(metadata$reads_type[index], pattern = type, ignore.case = TRUE)
    }))
    if (sum(is_seq_type) == 0) {
        stop(call. = FALSE, paste0(
            'The value in the "reads_type" column on row ', index, ' does not contain a known sequence type. ',
            'One of the following words must appear (case insensitive):\n',
            paste0('"', known_read_types, '"', collapse = ', '), '\n'
        ))
    }
    if (sum(is_seq_type) > 1) {
        stop(call. = FALSE, paste0(
            'The value in the "reads_type" column on row ', index, ' contains the names of multiple sequencing types. ',
            'Exactly one of the following words must appear (case insensitive):\n',
            paste0('"', known_read_types, '"', collapse = ', '), '\n'
        ))
    }
    return(known_read_types[is_seq_type])
}))

# Add user-supplied data as columns with modified names
cols_to_add <- colnames(metadata_original_samp)[colnames(metadata_original_samp) %in% known_columns_samp]

#since metadata_original_samp may have been modified, using its column numbers as index for original user input
unmodified_data <- unmodified_data[, colnames(metadata_original_samp) %in% known_columns_samp, drop= FALSE]

colnames(unmodified_data) <- paste0(user_column_name_prefix, colnames(unmodified_data))
metadata <- cbind(metadata, unmodified_data)

# Write output metadata
write.csv(metadata, file = args$output_path_samp, row.names = FALSE, na = '')
