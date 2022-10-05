#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

library(Seurat)
library(seuratTools)
library(infercnv)

normal_reference_mat <- readRDS("output/infercnv/reference_counts.rds")
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
	RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

make_seus_from_cellranger <- function(sample_path, normal_seu){
	# browser()
	mypath <- fs::path(sample_path)

	count_mat <- Seurat::Read10X(mypath)

	seu <- Seurat::CreateSeuratObject(count_mat, assay = "gene") %>%
		RenameCells(new.names = str_replace(colnames(.), "-", "."))

	seu_path <- path("output/seurat", paste0(path_file(sample_path), "_seu.rds"))

	seu <- seuratTools::clustering_workflow(seu, resolution = c(0.2, 0.4))

	saveRDS(seu, seu_path)

	# seu_merged <- merge(seu, normal_seu)
	#
	# seu_merged <- infercnv::add_to_seurat(seu_merged, fs::path("output/infercnv", path_file(sample_path)))
	#
	# seu_w_cnv <- seu_merged[,!grepl("normal", colnames(seu_merged))]
	#
	# seu_w_cnv <- seuratTools::clustering_workflow(seu_w_cnv, resolution = c(0.2, 0.4))
	#
	# seu_cnv_path <- path("output/seurat", paste0(path_file(sample_path), "_cnv_seu.rds"))
	#
	# saveRDS(seu_w_cnv, seu_cnv_path)
	#
	# seu_cnv_path

}

sample_paths <-
	fs::dir_ls("data/", glob = "*SRR*")

make_seus_from_cellranger(sample_paths[[1]], normal_seu)

purrr::map(sample_paths, make_seus_from_cellranger, normal_seu)

infercnv_obj <- readRDS("output/infercnv/SRR13884240/run.final.infercnv_obj")


seu_40 <- readRDS("output/seurat/SRR13884240_seu.rds") %>%
	RenameCells(new.names = str_replace(colnames(.), "-", "."))

normal_reference_mat <- readRDS("output/infercnv/reference_counts.rds")
normal_seu <- Seurat::CreateSeuratObject(normal_reference_mat) %>%
	RenameCells(new.names = paste0("normal_", str_replace(colnames(.), "-", ".")))

seu0 <- merge(seu_40, normal_seu)

seu_cells <- rownames(seu0@meta.data)

infercnv_cells <- colnames(infercnv_obj@expr.data)


# debug(infercnv::add_to_seurat)
seu1 <- infercnv::add_to_seurat(seu0, "output/infercnv/SRR13884240/")

seu2 <- seu1[,!grepl("normal", colnames(seu1))]

seu2 <- seurat_preprocess(seu2) %>%
	RunPCA() %>%
	RunUMAP(dims = 1:30)

has_cnv_cols <- str_subset(colnames(seu2@meta.data), "has_cnv*")

cnv_plots <- purrr::map(has_cnv_cols, ~DimPlot(seu2, group.by = .x))

pdf("~/tmp/cnvplots.pdf")
cnv_plots
dev.off()
