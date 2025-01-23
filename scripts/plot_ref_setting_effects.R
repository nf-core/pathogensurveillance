library(tibble)
library(purrr)
library(ggplot2)
library(dplyr)

raw_data <- system('wc -l path_surveil_data/assembly_metadata/*.tsv', intern = T)
raw_data <- raw_data[-length(raw_data)]

parsed_data <- tibble(
    line_count = as.numeric(map_chr(strsplit(raw_data, split = ' +'), `[`, 2)) - 1,
    file_path = map_chr(strsplit(raw_data, split = ' +'), `[`, 3),
    file_name = sub(basename(file_path), pattern = '.tsv$', replacement = ''),
    family = map_chr(strsplit(file_name, split = '--'), `[`, 1),
    setting = map_chr(strsplit(file_name, split = '--'), function(x) paste0(x[-1], collapse = '--')),
)
parsed_data$setting[parsed_data$setting == ''] <- 'default'

parsed_data %>%
    filter(line_count > 0) %>%
    ggplot(aes(x = setting, y = line_count, group = family, color = family)) +
    # geom_boxplot() +
    geom_line() +
    coord_trans(y = "log10")

