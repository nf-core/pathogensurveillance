# Deletes all task working directories in the "work" directory that had a non-zero exit code.


# Check that work directory exists
if (! dir.exists('work')) {
    stop(call. = FALSE, 'There is no "work"directory in the current working directory. Run this script from the location where a Nextflow pipeline was run.')
}

# Find all exit status files
exitcode_paths <- list.files("work", pattern = '.exitcode', recursive = TRUE, all.files = TRUE, full.names = TRUE)

# Read exit codes
exitcodes <- vapply(exitcode_paths, function(p) {
    out <- readLines(p)
    if (length(out) == 0) {
        return(NA_character_)
    } else {
        return(out)
    }
}, FUN.VALUE = character(1))

# Filter for non-zero exit codes
non_zero_exitcodes <- exitcodes[is.na(exitcodes) | exitcodes != "0"]

# Delete directories with non-zero exit codes
for (path in names(non_zero_exitcodes)) {
    cat(paste0('rm -rf ', dirname(path), '\n'))
}
