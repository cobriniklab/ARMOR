#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

library(Seurat)
library(seuratTools)
library(infercnv)

normal_reference_path <- "~/Homo_sapiens/infercnv/reference_mat.rds"

normal_reference_mat <- readRDS(normal_reference_path)
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
  RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

make_seus_from_cellranger <- function(sample_path){
  # browser()
  seu_path <- path("output/seurat", paste0(path_file(sample_path), "_seu.rds"))

  print(seu_path)

  mypath <- fs::path(sample_path, "outs/filtered_feature_bc_matrix")

  count_mat <- Seurat::Read10X(mypath)

  seu <- Seurat::CreateSeuratObject(count_mat, assay = "gene") %>%
    RenameCells(new.names = str_replace(colnames(.), "-", "."))

  saveRDS(seu, seu_path)
  print(glue::glue("saved {seu_path}"))

  return(seu)

}

append_infercnv_to_seu <- function(sample_id, normal_seu){
  # browser()
  seu_path <- path("output/seurat", paste0(path_file(sample_id), "_seu.rds"))

  seu <- readRDS(seu_path)

  seu_merged <- merge(seu, normal_seu)

  infercnv_dir <- fs::path("output/infercnv", path_file(sample_id))

  infercnv_obj <- readRDS(path(infercnv_dir, "/run.final.infercnv_obj"))

  seu_merged <- seu_merged[,colnames(seu_merged) %in% colnames(infercnv_obj@expr.data)]

  seu_merged <- infercnv::add_to_seurat(seu_merged, infercnv_dir)

  seu_w_cnv <- seu_merged[,!grepl("normal", colnames(seu_merged))]

  seu_cnv_path <- path("output/seurat", paste0(path_file(sample_id), "_cnv_seu.rds"))

  seu_w_cnv <- seuratTools::clustering_workflow(seu_w_cnv, resolution = c(0.2, 0.4))

  saveRDS(seu_w_cnv, seu_cnv_path)

  return(seu_w_cnv)

}

cellranger_paths <-
  fs::dir_ls("output/cellranger/", glob = "*SRR*") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seus <- purrr::map(cellranger_paths, make_seus_from_cellranger)

seu_paths <-
  cellranger_paths <-
  fs::dir_ls("output/seurat/", glob = "*SRR*") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seus <- purrr::map(seu_paths, readRDS)

marker_plots <- purrr::map(seus, ~plot_markers(.x, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE))

umap_plots <- purrr::imap(seus, ~(DimPlot(.x, group.by = "gene_snn_res.0.2") + labs(title = .y)))

library(patchwork)

patchworks <- purrr::map2(umap_plots, marker_plots, ~.x + .y + plot_layout(widths = c(3,1)))

pdf("results/patchworks.pdf", width = 10, height = 8)
patchworks
dev.off()

# exclude code ------------------------------

sample_ids <- c("SRR14800534", "SRR14800535", "SRR14800536", "SRR14800537", "SRR14800538", "SRR14800539", "SRR14800540", "SRR14800541", "SRR14800542", "SRR14800543")

test0 <- append_infercnv_to_seu("SRR17960480", normal_seu)

