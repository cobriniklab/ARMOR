#!/usr/bin Rscript

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

# Rprof(rprof_out)

library(numbat)
library(Seurat)
library(readr)
library(magrittr)
conflicted::conflict_prefer("rowSums", "Matrix")

nb = Numbat$new(out_dir = numbat_dir)

nb_path = paste0(numbat_dir, "_numbat.rds")

saveRDS(nb, nb_path)