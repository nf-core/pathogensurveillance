# If being run manually (no params), then get default params from yaml header
if (! exists("params")) {
    header <- rmarkdown::yaml_front_matter('index.qmd')
    params <- as.list(header$params)
}

# Parse metadata
group <- params$group
refs <- strsplit(params$refs, ';', fixed = TRUE)[[1]]
samp_data <- read_csv(params$samp_data, show_col_types = FALSE)
ref_data <- read_tsv(params$ref_data, col_types = 'dcccccccccccccccddc')

# Parse sendsketch data
sketch_data <- map_dfr(list.files(params$sendsketch), function(path) {
    data <- read_tsv(file.path(params$sendsketch, path), skip = 2,
                     show_col_types = FALSE)
    id <- sub(path, pattern = '\\.txt$', replacement = '')
    return(bind_cols(sample_id = rep(id, nrow(data)), data))
})

# Parse variant data
snp_tree_paths <- list.files(params$variant_data, pattern = "\\.treefile$", full.names = TRUE)
vcf_paths <- list.files(params$variant_data, pattern = "\\.vcf\\.gz$", full.names = TRUE)
snp_align_paths <- list.files(params$variant_data, pattern = "\\.fasta$", full.names = TRUE)

# Parse ANI matrix
ani_matrix <- read.csv(params$ani_matrix, check.names = FALSE)
rownames(ani_matrix) <- colnames(ani_matrix)

# Parse core gene phylogeny
core_phylo_path = params$core_phylo

# Parse quality control data
multiqc_report_path <- file.path(params$multiqc, 'multiqc_report.html')
quast_ref_names <- list.files(params$quast)
quast_report_paths <- file.path(quast_ref_names, 'report.html')

# Parse version data
raw_version_data <- unlist(read_yaml(params$versions))
version_data <- tibble(
    module = map_chr(strsplit(names(raw_version_data), split = '.', fixed = TRUE), `[`, 1),
    program = map_chr(strsplit(names(raw_version_data), split = '.', fixed = TRUE), `[`, 2),
    version = unname(raw_version_data)
)
