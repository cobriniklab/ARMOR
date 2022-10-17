#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(numbat)
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

numbat_seu_paths <-
  fs::dir_ls("output/seurat/", glob = "*SRR*numbat*") %>%
  purrr::set_names(str_extract(path_file(.), "SRR[0-9]*"))

seus <- purrr::map(numbat_seu_paths, readRDS)

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

make_phylo_heatmap <- function(numbat_dir){
    nb = Numbat$new(out_dir = numbat_dir)
    
    nb$plot_consensus()
    
    mypal = scales::hue_pal()(n_distinct(nb$clone_post$clone_opt))
    
    names(mypal) <- seq(1:length(mypal))
    
    numbat_phylo_heatmap <- nb$plot_phylo_heatmap(
        clone_bar = TRUE,
        p_min = 0.5,
        pal_clone = mypal
    )
    
    return(numbat_phylo_heatmap)
}

numbat_dirs = path("output/numbat", str_extract(path_file(numbat_seu_paths), "SRR[A-Z]*[0-9]*"))

phylo_heatmaps <- map(numbat_dirs, make_phylo_heatmap)

numbat_phylo_images <- dir_ls("output/numbat", glob = "*panel_2.png", recurse = TRUE) %>%
    purrr::map(read_image_as_plot) %>%
    identity()

infercnv_plots <- purrr::map(seus, FeaturePlot, features = cnv_cols)

numbat_cnv_plots <- purrr::map(seus, DimPlot, group.by = "clone_opt")

marker_plots <- purrr::map(seus, ~plot_markers(.x, metavar = "gene_snn_res.0.2", marker_method = "presto", return_plotly = FALSE))

umap_plots <- purrr::imap(seus, ~(DimPlot(.x, group.by = "gene_snn_res.0.2") + labs(title = .y)))

library(patchwork)

assemble_cnv_patchwork <- function(umap, markers, cnv_image, cnv_plot){
  mypatchwork = (umap + markers + cnv_image + cnv_plot) + plot_layout(ncol = 2, widths = c(2,2))

  return(mypatchwork)
}

patchworks <- purrr::pmap(list(umap_plots, marker_plots, infercnv_images, infercnv_plots), assemble_cnv_patchwork)


pdf("results/patchworks.pdf", width = 14, height = 14)
patchworks
dev.off()
