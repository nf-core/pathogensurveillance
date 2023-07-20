library(tidyverse)
library(ggplot2)
library(plyr)
library(dplyr)

###### Load FILES

### PIRATE.gene_families.ordered.originalIDs.tsv
gene_fam <- read.csv(file = paste0(path,"PIRATE.gene_families.ordered.originalIDs.tsv"), 
                     sep = "\t", stringsAsFactors = FALSE, fill = TRUE, header = TRUE)


### CORE genes

## Subset core genes with only one copy (no paralogs), from gene_fam table (original input)
core_complete_tab <- gene_fam %>% 
  filter(number_genomes == core & max_dose ==1)

for (r in 1:ncol(core_complete_tab)) {
  tmp_rm <- str_remove_all(core_complete_tab[,r], pattern = "cds-") # remove pattern "cds-"
  #tmp_rm <- str_remove_all(tmp_rm, pattern = ",")   # remove pattern ","
  core_complete_tab[,r] <- as.matrix(tmp_rm)
}

## Save subset table to file
write.table(core_complete_tab, file = "core_complete_tab.txt", col.names = TRUE, row.names = F, sep= "\t", quote = FALSE)


# Get gene families of core genes
gf_names <- core_complete_tab$gene_family
write.table(gf_names, file = "core_geneFam_list.txt", col.names = F, row.names = F, sep= "\t", quote = FALSE)

