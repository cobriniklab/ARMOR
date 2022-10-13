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

cellranger_paths <-
  fs::dir_ls("output/cellranger/", glob = "*SRR*") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seu_cnv_paths <-
  fs::dir_ls("output/seurat/", glob = "*SRR*cnv*") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seus <- purrr::map(seu_cnv_paths, readRDS)

cnv_cols <- c("proportion_dupli_1", "proportion_dupli_2", "proportion_dupli_6", "proportion_loss_16")

read_image_as_plot <- function(image_path){
  infercnv_image <- magick::image_read(image_path)

  infercnv_image <- ggplot() +
    ggpubr::background_image(infercnv_image) +
    # coord_fixed() +
    NULL

  return(infercnv_image)
}

infercnv_image_paths <- dir_ls("output/infercnv/", glob = "*infercnv.png", recurse = TRUE)

dir_create("results/infercnv")

new_infercnv_image_paths <- path("results/infercnv", paste0(path_file(path_dir(infercnv_image_paths)), "_infercnv.png"))

map2(infercnv_image_paths, new_infercnv_image_paths, fs::file_copy)

infercnv_images <- dir_ls("results/infercnv", glob = "*.png") %>%
  purrr::map(read_image_as_plot) %>%
  identity()

cnv_plots <- purrr::map(seus, FeaturePlot, features = cnv_cols)

marker_plots <- purrr::map(seus, ~plot_markers(.x, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE))

umap_plots <- purrr::imap(seus, ~(DimPlot(.x, group.by = "gene_snn_res.0.2") + labs(title = .y)))

library(patchwork)

assemble_cnv_patchwork <- function(umap, markers, cnv_image, cnv_plot){
  mypatchwork = (umap + markers + cnv_image + cnv_plot) + plot_layout(ncol = 2, widths = c(2,2))

  return(mypatchwork)
}

patchworks <- purrr::pmap(list(umap_plots, marker_plots, infercnv_images, cnv_plots), assemble_cnv_patchwork)


pdf("results/patchworks.pdf", width = 14, height = 14)
patchworks
dev.off()

# exclude code ------------------------------

test0 <- append_infercnv_to_seu("SRR17960480", normal_seu)

