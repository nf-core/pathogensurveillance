#!/usr/bin/env Rscript

# MIT License
#
# Copyright (c) Zachary S.L. Foster and Niklaus J. Grunwald
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
