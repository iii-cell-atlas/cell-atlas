#library(SeuratData)
library(shiny)
library(shinythemes)
library(plotly)
library(Seurat)
#install.packages("magrittr")
library(magrittr)
library(dplyr)
library(cowplot)
#install.packages("ggplot2")
library("ggplot2")
library(rlang)
library(tidyr)
library(tibble)
library(tidyverse)
library(ggrepel)
library(grid)
library(data.table)
library(shinycssloaders)
library(DT)

# adapt to your path
#setwd("D:/GCRF_UoG/Vicky_JCR_Shiny/Tcells_amp_adopted/")

##Reading in the list of precomputed Seurat objects (Resolution 0.15 0.25&0.35, 0.45&0.55 computed separately read in and combined then written to disk and )
tcells_combined_umap_list_res<-readRDS("tcells_combined_umap_list_res.rds")
# tcells_combined_clusters_tables_res1<-readRDS("D:/GCRF_UoG/Vicky_JCR_Shiny/tcells_combined_clusters_tables_res.rds")
# saveRDS(tcells_combined_clusters_tables_res1, "D:/GCRF_UoG/Vicky_JCR_Shiny/tcells_combined_clusters_tables_res1.rds")
# tcells_combined_clusters_tables_res2_3<-readRDS("D:/GCRF_UoG/Vicky_JCR_Shiny/tcells_combined_clusters_tables_res2_3.rds")
# tcells_combined_clusters_tables_res4_5<-readRDS("D:/GCRF_UoG/Vicky_JCR_Shiny/tcells_combined_clusters_tables_res4_5.rds")
# tcells_combined_clusters_tables_res = c(tcells_combined_clusters_tables_res1, tcells_combined_clusters_tables_res2_3, tcells_combined_clusters_tables_res4_5)
# saveRDS(tcells_combined_clusters_tables_res, "D:/GCRF_UoG/Vicky_JCR_Shiny/tcells_combined_clusters_tables_res.rds")

##Reading in the list of precomputed table of cluster markers
tcells_combined_clusters_tables_res<-readRDS("tcells_combined_clusters_tables_res.rds")

##Reading in the list of precomputed table of differential expressed (DE) genes
tcells_combined_de_tables = readRDS("tcells_combined_de_tables.rds")

##Reading in the list of precomputed table of cluster markers
tcells_combined_de_ggplots_table = readRDS("tcells_combined_de_ggplots_table.rds")

##Reading in all the genes that are present in WT1, WT2 and KO raw objects
all_genes_common_in_all_groups = readRDS("all_genes_common_in_all_groups.rds")

##Reading in and processing table of uniprot links for all mouse genes
# uniprot_info_raw = fread("/GCRF_UoG/Vicky_JCR_Shiny/uniprot table/unipro-mouseID")
# uniprot_info_raw$uniprot = paste('<a href="https://www.uniprot.org/uniprot/',uniprot_info_raw$Entry,'" target="_blank">', uniprot_info_raw$Entry,'</a>', sep = "")
# write.table(uniprot_info_raw, "/GCRF_UoG/Vicky_JCR_Shiny/uniprot_info_with_link", row.names = F, sep = "\t", quote = F)
uniprot_info = fread("uniprot_info_with_link", stringsAsFactors = F)

##Declaring and assigning variables
dim=15

res1 = 0.15
res2 = 0.55
diff_res = 0.1
cluster_names = c("Il17a +ve cells", "Ccr7 +ve cells", "Ly6c2 +ve cells", "Gzma +ve cells", "Cdk6 +ve cells")
fav_genes = c("Cdk6", "Gzma", "Ly6c2", "Ccr7", "Il17a")
conditions = c("CD18 KO", "WT")

##Code block for processing raw 10X data and computing the Seurat objects and cluster marker tables
#if (!(all(exists("tcells.combined.umap.list"), exists("tcells.combined.clusters.tables")))) {
if (!(all(exists("tcells_combined_umap_list_res"), exists("tcells_combined_clusters_tables_res"), exists("tcells_combined_de_tables")))) {
  #Merge first the two WT 1 and 3
  WT1.data <- Read10X(data.dir = "WT1")
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
  KO.data <- Read10X(data.dir = "KO1")
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
  
  WT3.data <- Read10X(data.dir = "WT3")
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
  
  
  ##Precomputing and list of seurat objects with different clusters allowing change of PCs
  # tcells.combined.umap.list = lapply(seq(dim1, dim2, by = diff_dim), function(x) {
  #   a = RunUMAP(Combined.filt, reduction = "pca", dims = 1:x)
  #   lapply(seq(dim1, dim2, by = diff_dim), function(x) {
  #     b = FindNeighbors(a, reduction = "pca", dims = 1:x)
  #     lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(b, resolution = x))
  #   })
  # })
  
  ##Precomputing and saving the list of Seurat objects with different clusters through adjusting of resolution from 0.15, 0.25, 0.35, 0.45 & 0.55
  tcells_combined_umap_list_res = lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(Combined.filt, resolution = x))
  saveRDS(tcells_combined_umap_list_res, "tcells_combined_umap_list_res.rds")
  
  ##Precomputing and saving the list of Seurat objects with different clusters through adjusting of resolution from 0.15, 0.25, 0.35, 0.45 & 0.55
  # tcells.combined.clusters.tables = lapply(tcells.combined.umap.list, function(x) { 
  #   lapply(x, function(x) {
  #     lapply(x, function(x){
  #       DefaultAssay(x) = "RNA"
  #       #cluster.markers = lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
  #       cluster.markers = lapply(0:1, function(y) {
  #         FindConservedMarkers(x, ident.1 = y, grouping.var = "sample")
  #         
  #       })
  #     }) 
  #   }) 
  # })
  # 
  # saveRDS(tcells.combined.clusters.tables, "tcells.combined.clusters.tables.rds")
  
  ##Precomputing and saving the list of tables of cluster markers
  tcells_combined_clusters_tables_res = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      FindConservedMarkers(x, ident.1 = y, grouping.var = "group")
      
    })
  })
  
  saveRDS(tcells_combined_clusters_tables_res, "tcells_combined_clusters_tables_res.rds")
  
  
  ##Precomputing and saving the list of tables of DE genes per cluster
  tcells_combined_de_tables = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    x$celltype.group <- paste(Idents(x), x$group, sep = "_")
    x$celltype <- Idents(x)
    Idents(x) <- "celltype.group"
    
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      FindMarkers(x, ident.1 = paste(y, "KO", sep = "_"), ident.2 = paste(y, "WT", sep = "_"), verbose = T, min.cells.group = 2)
      
    })
  })
  
  saveRDS(tcells_combined_de_tables, "tcells_combined_de_tables.rds")
  
  ##Precomputing and saving the list of ggplot per cluster for all resolutions
  tcells_combined_de_ggplots_table = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      cells_type <- subset(x, idents = y)
      #Idents(cells_type) <- "sample"
      Idents(cells_type) <- "group"
      avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
      avg.cells$gene <- rownames(avg.cells)
      avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene)) %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = "Gene names  (primary )")) %>% dplyr::distinct(., gene, .keep_all = T)%>% select(gene, KO, WT , `Protein names`, uniprot)
      
    })
  })
  
  saveRDS(tcells_combined_de_ggplots_table, "tcells_combined_de_ggplots_table.rds")
  
  # ##Precomputing and saving the list of ggplot per cluster for all resolutions
  # tcells_combined_de_ggplots = lapply(tcells_combined_de_ggplots_table, function(x) { 
  #   DefaultAssay(x) = "RNA"
  #   
  #   lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
  #     theme_set(theme_cowplot())
  #     grob <- grobTree(textGrob("Click on points to diplay more information about the gene", x=0.1,  y=0.95, hjust=0,
  #                               gp=gpar(col="red", fontsize=9, fontface="italic")))
  #     ggplot(data=avg.cells, aes_string("KO","WT")) + geom_point() + annotation_custom(grob)
  #     
  #   })
  # })
  # 
  # saveRDS(tcells_combined_de_tables, "tcells_combined_de_ggplot.rds")
  
}

##functions


which_numeric_cols = function(dat) {
  which(sapply(seq(ncol(dat)), function(i) {
    is.numeric(dat[,i])
  }))
}

