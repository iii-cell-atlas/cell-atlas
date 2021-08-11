#!/usr/bin/env Rscript
library(dplyr)
library(httr)
library(MAST)
library(unixtools)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
library(ggmin)
library(pheatmap)
library(ggrepel)
library(shiny)
library(shinythemes)
library(plotly)
library(magrittr)
library(rlang)
library(tidyr)
library(tibble)
library(tidyverse)
library(grid)
library(data.table)
library(shinycssloaders)
library(DT)
library(shinydashboard)
library(shinyjs)

setwd(paste0(getwd(), "/data/"))

if (all(
  file.exists("combined_umap_list_res_skinny.rds"),
  file.exists("combined_clusters_tables_res.rds"),
  file.exists("combined_de_tables.rds"),
  file.exists("combined_de_ggplots_table.rds"),
  file.exists("all_genes_common_in_all_groups.rds"))) {
  
  ##Reading in the list of precomputed Seurat objects (Resolution 0.15 0.25&0.35, 0.45&0.55 computed separately read in and combined then written to disk and )
  combined_umap_list_res_skinny <- readRDS("combined_umap_list_res_skinny.rds")
  
  ##Reading in the list of precomputed table of cluster markers
  combined_clusters_tables_res <- readRDS("combined_clusters_tables_res.rds")
  
  ##Reading in the list of precomputed table of differential expressed (DE) genes
  combined_de_tables <- readRDS("combined_de_tables.rds")
  
  ##Reading in the list of precomputed table of cluster markers
  combined_de_ggplots_table <- readRDS("combined_de_ggplots_table.rds")
  
  ##Reading in all the genes that are present in WT1, WT2 and KO raw objects
  all_genes_common_in_all_groups <- readRDS("all_genes_common_in_all_groups.rds")
  
}

#~~SCRIPT-SPECIFIC~~#
##Reading in and processing table of uniprot links
uniprot_info = data.frame()


#~~SCRIPT-SPECIFIC~~#
##Declaring and assigning variables
#vvv Used in this script only:
dim <- 10
#^^^
#vvv Used in this script and app.R:
res1 <- 0.1
res_default <- 0.4
res2 <- 1.0
diff_res <- 0.1
#^^^
#vvv Used in app.R only:
cluster_names <- c("")
fav_genes <- c("")
conditions <- c("")
umap_names <- c("")
cluster.colours <- c("")
group.cols <- c("")

names(group.cols) <- conditions
# Group up conditions into combinations of 2
pairwise <- combn(conditions, 2)
conds = lapply(1:ncol(pairwise), function(x) paste(pairwise[,x], collapse = " VS ")) %>% unlist()
#^^^

#vvv Only the following remaining objects are required by app.R:
# tcells_combined_umap_list_res_skinny
# modify_stop_propagation (a function)
# tcells_combined_de_ggplots_table
# tcells_combined_de_tables
# tcells_combined_clusters_tables_res
# all_genes_common_in_all_groups
#^^^

if (!(all(
  exists("combined_umap_list_res_skinny"),
  exists("combined_clusters_tables_res"),
  exists("combined_de_tables")))) {
  
  #Merge first the two WT 1 and 3
  WT1.data <- Read10X(data.dir = "t_cell_raw_data/WT1")
  WT1 <- CreateSeuratObject(counts = WT1.data, project = "WT1", min.cells = 3)
  WT1$sample <- "WT1"
  WT1$group <- "WT"
  WT1[["percent.mt"]] <- PercentageFeatureSet(object = WT1, pattern = "^mt-")
  plot1 <- FeatureScatter(object = WT1, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = WT1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  WT1 <- subset(WT1, subset = nFeature_RNA > 500 & nFeature_RNA < 2200 & percent.mt < 18)
  
  # store mitochondrial percentage in object meta data
  WT1 <- PercentageFeatureSet(WT1, pattern = "^mt-", col.name = "percent.mt")
  
  
  
  # run sctransform
  WT1 <- SCTransform(WT1, vars.to.regress = "percent.mt", verbose = FALSE)
  
  
  ##capturing all genes at this stage
  all_genes_wt1 = rownames(WT1)
  saveRDS(all_genes_wt1, "all_genes_wt1.rds")
  ##capturing all genes at this stage
  
  ### load the KO data
  KO.data <- Read10X(data.dir = "t_cell_raw_data/KO1")
  KO <- CreateSeuratObject(counts = KO.data, project = "KO1", min.cells = 3)
  KO$sample <- "KO1"
  KO$group <- "KO"
  KO[["percent.mt"]] <- PercentageFeatureSet(object = KO, pattern = "^mt-")
  plot1 <- FeatureScatter(object = KO, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = KO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  # store mitochondrial percentage in object meta data
  KO <- PercentageFeatureSet(KO, pattern = "^mt-", col.name = "percent.mt")
  
  # run sctransform
  KO <- SCTransform(KO, vars.to.regress = "percent.mt", verbose = FALSE)
  
  ##capturing all genes at this stage
  all_genes_ko = rownames(KO)
  saveRDS(all_genes_ko, "all_genes_ko.rds")
  ##capturing all genes at this stage
  
  WT3.data <- Read10X(data.dir = "t_cell_raw_data/WT3")
  WT3 <- CreateSeuratObject(counts = WT3.data, project = "WT3", min.cells = 3)
  WT3$sample <- "WT2"
  WT3$group <- "WT"
  WT3[["percent.mt"]] <- PercentageFeatureSet(object = WT3, pattern = "^mt-")
  plot1 <- FeatureScatter(object = WT3, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(object = WT3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  WT3 <- subset(WT3, subset = nFeature_RNA > 500 & nFeature_RNA < 3200 & percent.mt < 18)
  
  
  # store mitochondrial percentage in object meta data
  WT3 <- PercentageFeatureSet(WT3, pattern = "^mt-", col.name = "percent.mt")
  
  # run sctransform
  WT3 <- SCTransform(WT3, vars.to.regress = "percent.mt", verbose = FALSE)
  
  ##capturing all genes at this stage
  all_genes_wt3 = rownames(WT3)
  saveRDS(all_genes_wt3, "all_genes_wt3.rds")
  ##capturing all genes at this stage
  
  WT3.downsample = subset(WT3, cells = sample(Cells(WT3), 4000))
  
  ##capturing all genes at this stage
  all_genes_wt3_ds = rownames(WT3)
  saveRDS(all_genes_wt3_ds, "all_genes_wt3_ds.rds")
  ##capturing all genes at this stage
  
  
  
  ## estimated in our script the best value here
  #WT.combined <- JackStraw(object = WT.combined, num.replicate = 100, dims=30)
  #WT.combined <- ScoreJackStraw(object = WT.combined, dims = 1:30)
  #JackStrawPlot(object = WT.combined, dims = 1:30)
  #ElbowPlot(object = WT.combined, ndims = 40)
  # merge the two WT's
  #dim=15
  
  TC.anchors <- FindIntegrationAnchors(object.list = list(WT1,KO), dims = 1:dim)
  Combined <- IntegrateData(anchorset = TC.anchors, dims = 1:dim)
  DefaultAssay(Combined) <- "integrated"
  Combined <- ScaleData(Combined, verbose = FALSE)
  Combined <- RunPCA(Combined, npcs = dim, verbose = FALSE)
  Combined <- RunUMAP(Combined, reduction = "pca", dims = 1:dim)
  Combined <- FindNeighbors(Combined, reduction = "pca", dims = 1:dim)
  Combined <- FindClusters(Combined, resolution = 0.2)
  p1 <- DimPlot(Combined, reduction = "umap", group.by = "sample")
  p2 <- DimPlot(Combined, reduction = "umap", label = TRUE)
  
  #x11()
  plot_grid(p1, p2)
  Combined[["UMI"]] <-  Combined$nCount_RNA  # Why divided by 100
  Combined[["genes"]] <-  Combined$nFeature_RNA
  FeaturePlot(Combined, features = "UMI")
  
  ##Number of PCs selected based on prior expert analysis
  dim=15
  WT.anchors <- FindIntegrationAnchors(object.list = list(Combined, WT3.downsample), dims = 1:dim)
  Three.combined <- IntegrateData(anchorset = WT.anchors, dims = 1:dim)
  DefaultAssay(Three.combined) <- "integrated"
  Three.combined <- ScaleData(Three.combined, verbose = FALSE)
  Three.combined <- RunPCA(Three.combined, npcs = dim, verbose = FALSE)
  Three.combined <- RunUMAP(Three.combined, reduction = "pca", dims = 1:dim)
  Three.combined <- FindNeighbors(Three.combined, reduction = "pca", dims = 1:dim)
  Three.combined <- FindClusters(Three.combined, resolution = 0.1)
  p1 <- DimPlot(Three.combined, reduction = "umap", group.by = "sample")
  p2 <- DimPlot(Three.combined, reduction = "umap", label = TRUE, label.size = 5)
  plot_grid(p1, p2)
  
  # filtering away other clusters
  Combined.filt <- subset(Three.combined, idents = c("0","1","2"), invert = FALSE)
  DimPlot(Combined.filt, reduction = "umap", split.by = "sample")
  
  # re-clusterting, without the other cell types
  DefaultAssay(object = Combined.filt) <- "integrated"
  # Run the standard workflow for visualization and clustering
  Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
  Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
  
  # t-SNE and Clustering
  Combined.filt <- RunUMAP(object = Combined.filt, reduction = "pca", dims = 1:dim)
  Combined.filt <- FindNeighbors(object = Combined.filt, reduction = "pca", dims = 1:dim)
  #Combined.filt <- FindClusters(Combined.filt, resolution = 0.15)
  
  
  ### decrease the amount of clusters, and split the other one (green, killers)
  # filtering away other clusters
  # Combined.filt2 <- subset(Combined.filt, idents = c("0","1","2","3","5"), invert = FALSE)
  # 
  # Combined.filt<-Combined.filt2
  # # re-clusterting, without the other cell types
  # DefaultAssay(object = Combined.filt) <- "integrated"
  # # Run the standard workflow for visualization and clustering
  # Combined.filt <- ScaleData(object = Combined.filt, verbose = FALSE)
  # Combined.filt <- RunPCA(object = Combined.filt, npcs = dim, verbose = FALSE)
  # 
  all_genes_common_in_all_groups = Reduce(intersect,list(all_genes_ko,all_genes_wt1,all_genes_wt3))
  saveRDS(all_genes_common_in_all_groups, "all_genes_common_in_all_groups.rds")
}
