#!/usr/bin/env Rscript
library(dplyr)
library(httr)
library(MAST)
library(unixtools)
library(Seurat)
library(ggplot2)
library(cowplot)
library(future)
#library(ggmin)
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


# adapt to your path
#setwd("/home/jr345y/tcell_cd18KO_compact/tcell_cd18KO_compact/")

##Reading in the list of precomputed Seurat objects (Resolution 0.15 0.25&0.35, 0.45&0.55 computed separately read in and combined then written to disk and )
tcells_combined_umap_list_res_skinny<-readRDS("tcells_combined_umap_list_res_skinny.rds")

##Reading in the list of precomputed table of cluster markers
tcells_combined_clusters_tables_res<-readRDS("tcells_combined_clusters_tables_res.rds")

##Reading in the list of precomputed table of differential expressed (DE) genes
tcells_combined_de_tables = readRDS("tcells_combined_de_tables.rds")

##Reading in the list of precomputed table of cluster markers
tcells_combined_de_ggplots_table = readRDS("tcells_combined_de_ggplots_table.rds")

##Reading in all the genes that are present in WT1, WT2 and KO raw objects
all_genes_common_in_all_groups = readRDS("all_genes_common_in_all_groups.rds")

##Reading in and processing table of links to PlasmoDB.. 
uniprot_info = fread("Pberghei.Anntation.txt", stringsAsFactors = F)


##Declaring and assigning variables
dim=24

res1 = 0.51
res2 = 0.51
diff_res = 0.2
#cluster_names = c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
#                  "Early Trophs 2", "Mid Trophs", "Late Trophs","Schizonts","Early Gams","Males","Females");
cluster_names = c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
                  "Early Trophs 2", "Mid Trophs", "Late Trophs","Schizonts","Early Gams","Males","Females");
#
#Early Trophs 2 Mid Rings      Mid Trophs     Early Rings    Early Trophs 1 Late Rings     Late Trophs    Early Gams    Females        Males          Schizonts 
fav_genes = c("PBANKA-1437500.1","PBANKA-1334300.1","PBANKA-0515000.1")
conditions = c("B", "BM","S")
umap_names = c("Early Rings", "Mid Rings", "Late Rings", "Early Trophs 1", 
                "Early Trophs 2", "Mid Trophs", "Late Trophs","Schizonts","Early Gams","Males","Females")

pairwise <- combn(conditions, 2)
conds = lapply(1:ncol(pairwise), function(x) paste(pairwise[,x], collapse = " VS ")) %>% unlist()
cluster.colours <-c("#A42537","dodgerblue2","#95E949","#FF1E1A")
group.cols <- c("red", "deepskyblue","green","blue")
names(group.cols) <- conditions
choice_gene = "PBANKA-1437500.1"
cond = "Mouse tissues"





if (!(all(exists("tcells_combined_umap_list_res_skinny"), exists("tcells_combined_clusters_tables_res"), exists("tcells_combined_de_tables")))) {
  

 
  
  Combined.filt<-readRDS("Pb.combined.noMCA.rds")
  DimPlot(Combined.filt)
  ## transfer the correct variables
  Combined.filt$sample<-Combined.filt$orig.ident; # orig.ident is the default name for a sample in Seurat
  
  ## in this case the group is the organ which is S for spleen BM for Bone marrow and B for blood
  Combined.filt$group<-Combined.filt$Organ
  Combined.filt$groups<-Combined.filt$Organ
  
  ## Potential bug, as some genes are not in all objects!
  all_genes_common_in_all_groups = rownames(Combined.filt@assays$RNA); #Reduce(intersect,list(all_genes_ko,all_genes_wt1,all_genes_wt3))
  saveRDS(all_genes_common_in_all_groups, "all_genes_common_in_all_groups.rds")

  Combined.filt <- FindNeighbors(Combined.filt)
  
  Combined.filt <-FindClusters(Combined.filt)
  ##Precomputing and saving the list of Seurat objects with different clusters through adjusting of resolution from 0.15, 0.25, 0.35, 0.45 & 0.55
  tcells_combined_umap_list_res <- lapply(seq(res1, res2, by = diff_res), function(x) FindClusters(Combined.filt, resolution = x))
  tcells_combined_umap_list_res[[2]] <- Combined.filt;
  tcells_combined_umap_list_res <-readRDS("tcells_combined_umap_list_res.rds")
  remove(tcells_combined_umap_list_res[[2]])
  tcells_combined_umap_list_res[[2]]=NULL
  
  ##Precomputing and saving the conserved markers
  tcells_combined_clusters_tables_res = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
      FindConservedMarkers(x, ident.1 = y, grouping.var = "group")
      
    })
  })
  
  saveRDS(tcells_combined_clusters_tables_res, "tcells_combined_clusters_tables_res.rds")
  
  
  ##Generating pairwise list for all DE group comparisons per cluster
  grps = unique(tcells_combined_umap_list_res[[1]]@meta.data$group)
  pairwise <- combn(grps, 2)
  
  ##Precomputing and saving the list of tables of DE genes per cluster
  tcells_combined_de_tables = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    x$celltype.group <- paste(Idents(x), x$group, sep = "_")
    x$celltype <- Idents(x)
    Idents(x) <- "celltype.group"
    
    lapply(unique(x$Stages), function(y) {
      lapply(1:ncol(pairwise), function(z) {
   #     print(y)
  #      print(z)
  #      print (paste(y, pairwise[1,z], sep = "_"))
        tryCatch(FindMarkers(x, ident.1 = paste(y, pairwise[1,z], sep = "_"), ident.2 = paste(y, pairwise[2,z], sep = "_"), 
                    verbose = T, min.cells.group = 3, assay = "RNA"), error=function(e) NULL)
        
      })
    })
  })
  
  saveRDS(tcells_combined_de_tables, "tcells_combined_de_tables.rds")
  
  ##Precomputing and saving the list of ggplot per cluster for all resolutions
  tcells_combined_de_ggplots_table = lapply(tcells_combined_umap_list_res, function(x) { 
    DefaultAssay(x) = "RNA"
    
#    lapply(0:(length(unique(x$seurat_clusters))-1), function(y) {
  # the seurat_cluster are just niumbers but we need the cluster names... 
        lapply(unique(x$Stages), function(y) {
      cells_type <- subset(x, idents = y)
      #Idents(cells_type) <- "sample"
      Idents(cells_type) <- "group"
      avg.cells <- log1p(AverageExpression(cells_type, verbose = FALSE)$RNA)
      avg.cells$gene <- rownames(avg.cells)
      avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene)) %>% dplyr::left_join(x = ., y = uniprot_info, by = c("gene" = "Gene")) %>% dplyr::distinct(., gene, .keep_all = T)%>% select(gene, `Protein name`, PlasmoDB, B, BM, S)
   #   avg.cells <- avg.cells %>% filter(!grepl("^mt-", gene))
    })
  })
  
  saveRDS(tcells_combined_de_ggplots_table, "tcells_combined_de_ggplots_table.rds")
  
  saveRDS(tcells_combined_umap_list_res,"tcells_combined_umap_list_res.rds")
  
  
  tcells_combined_umap_list_res_skinny = tcells_combined_umap_list_res
  tcells_combined_umap_list_res_skinny[[1]]@assays$integrated = NULL
  tcells_combined_umap_list_res_skinny[[1]]@assays$SCT = NULL
  tcells_combined_umap_list_res_skinny[[2]] =NULL
  saveRDS(tcells_combined_umap_list_res_skinny, "tcells_combined_umap_list_res_skinny.rds")
  
}

##functions
##function to make sidebar menu expanded my default
modify_stop_propagation <- function(x) {
  x$children[[1]]$attribs$onclick = "event.stopPropagation()"
  x
}
