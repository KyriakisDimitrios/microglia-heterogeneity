
source("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/Snakemake/2.Functions.R")

library(scater)
library(sctransform)
library(Seurat)
library(scran)
library(tictoc)
library(crayon)
library(cluster)
library(Routliers)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(jcolors)
library(NMF)
library(jcolors)
library(ggpubr)
library(cowplot)

set.seed(123)
tool="seurat"
project ="MBSYN"
dataset <- project
Data_select <- data_selection(project)

WORKDIR <- Data_select$WORKDIR
list_of_files <- Data_select$list_of_files
condition_names <- Data_select$condition_names
organism<- Data_select$organism
file<- Data_select$file
data_10x<- Data_select$data_10x
setwd(Data_select$WORKDIR)


color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_cells <- c(rgb(255,0,0, maxColorValue = 255),rgb(0,0,255, maxColorValue = 255),rgb(0,0,160, maxColorValue = 255),rgb(255,100,177, maxColorValue = 255),rgb(128,0,128, maxColorValue = 255),rgb(128,64,0, maxColorValue = 255),rgb(0,255,0, maxColorValue = 255))
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)


# ================================ SETTING UP ======================================== #
# Number of cells to use
imputation = FALSE
remove_mt=TRUE
remove_ribsomal=TRUE
n_cores=5
elbow = FALSE


setwd("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/")
NewDir <- paste0("Rescue_MBSYN_",tool,"_elbow_",elbow,"_remove_ribsomal_",remove_ribsomal)
dir.create(NewDir)
setwd(NewDir)





Combined <- readRDS("1.Combined.rds")
# ==============================================================================================
# ======================================= Clustering  ==========================================
# ==============================================================================================
dir.create("2.Clustering")
setwd("2.Clustering")
Combined <- reduce_dim(Combined,project=project)$Combined

plot_nFeatures(Combined,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)
plot_tot_mRNA(Combined,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)

pdf("QC.pdf")
FeaturePlot(Combined,features = "nFeature_RNA",reduction = "tsne",order=T)
FeaturePlot(Combined,features = "nCount_RNA",reduction = "tsne",order=T)
FeaturePlot(Combined,features = "nFeature_RNA",order=T)
FeaturePlot(Combined,features = "nCount_RNA",order=T)
dev.off()

pdf(paste("Umap.pdf",sep="_"))
DimPlot(object = Combined, reduction = "umap", group.by = "condition",cols = color_cond)
DimPlot(object = Combined, reduction = "umap", group.by = "Cluster", label = TRUE,cols = color_clust)
dev.off()

pdf(paste("TSNE.pdf",sep="_"))
DimPlot(object = Combined, reduction = "tsne", group.by = "condition",cols = color_cond)
DimPlot(object = Combined, reduction = "tsne", group.by = "Cluster", label = TRUE,cols = color_clust)
dev.off()

# ------------------------------------------------------------------------------------------------

setwd("../")

saveRDS(Combined,"2.Combined.rds")
