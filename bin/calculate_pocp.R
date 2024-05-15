#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
args <- as.list(args)
names(args) <- c('input_path', 'output_path')

piratetab <- read.table(args$input_path, header = TRUE, stringsAsFactors = FALSE, check.names=FALSE)

row.names(piratetab) <- piratetab$Gene
piratetab$Gene <- NULL
piratetab_tr <- as.data.frame(t(piratetab))
piratetab_tr$strain <- row.names(piratetab_tr)

pocp <- function(a, b) {
	  round(2*sum(a == '1' & b == '1')/(sum(a == '1') + sum(b == '1')) * 100, digits = 2)
}

tmp <- asplit(piratetab_tr, 1)
result <- outer(tmp, tmp, Vectorize(pocp))
dimnames(result) <- list(rownames(piratetab_tr), rownames(piratetab_tr))

write.table(result, file = args$output_path, sep = "\t", quote = FALSE, row.names = FALSE)
