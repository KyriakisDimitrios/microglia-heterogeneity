
source("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/Snakemake/3.Functions.R")

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





Combined <- readRDS("2.Combined.rds")





# =================================== CELL TYPES =================================================
dir.create("3.Cell_Types")
setwd("3.Cell_Types")
gene_list_name <- "Cell_Type"
file <- paste0(WORKDIR,"/Gene_Lists/Paper_Cell_types.txt")
res_cell_type_assignment <- cell_type_assignment(object=Combined,tab_name="Cell_Type",file=file)
Combined = res_cell_type_assignment$object
r_annot = res_cell_type_assignment$r_annot


# ============ ADVANCED QUALITY FILTERING =========================
cell_types <- Combined$Cell_Type
save_cell_types<-Combined$Cell_Type
cell_types[ which(Combined$Cluster %in% c(1))] <- "Astrocytes"
cell_types[ which(Combined$Cluster %in% c(2))] <- "Microglia"
cell_types[ which(Combined$Cluster %in% c(3))] <- "Oligodendrocytes"
cell_types[ which(Combined$Cluster %in% c(4))] <- "Endothelial"
cell_types[ which(Combined$Cluster %in% c(6,7))] <- "Ependymal"
cell_types[ which(Combined$Cluster %in% c(5,8))] <- "Neurons"
cell_types[ which(Combined$Cluster %in% c(9))] <- "Hybrid"
Combined$Cell_Type <-cell_types



gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
plot_bar_cells(gene_list,object=Combined,save=TRUE)


pie_per_CellType(Combined,target="Cell_Type")



# cell_types <- Combined$Cell_Type
# cell_types <- as.vector(Combined$Cluster)
# save_cell_types<-Combined$Cell_Type
# cell_types[ which(Combined$Cluster %in% c(1))] <- "Astrocytes"
# cell_types[ which(Combined$Cluster %in% c(2))] <- "Microglia"
# cell_types[ which(Combined$Cluster %in% c(3))] <- "Oligodendrocytes"
# cell_types[ which(Combined$Cluster %in% c(4))] <- "Endothelial"
# cell_types[ which(Combined$Cluster %in% c(5))] <- "Hybrid 1"
# cell_types[ which(Combined$Cluster %in% c(6))] <- "Ependymal"
# cell_types[ which(Combined$Cluster %in% c(7))] <- "Epithelial"
# cell_types[ which(Combined$Cluster %in% c(8))] <- "Neurons"
# cell_types[ which(Combined$Cluster %in% c(9))] <- "Hybrid 2"
# Combined$Cell_Type <-cell_types

color_clust_ord_cl <- color_clust[c(1,4,6,7,5,9,2,8,3)]
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_clust_ord_cl,State=color_clust)

# ================================================================

pdf(paste("Umap.pdf",sep="_"))
DimPlot(object = Combined, reduction = "umap", group.by = "condition",cols = color_cond)
DimPlot(object = Combined, reduction = "umap", group.by = "Cluster", label = TRUE,cols = color_clust)
DimPlot(object = Combined, reduction = "umap", group.by = "Cell_Type", label = TRUE,cols = color_cells)
dev.off()

pdf(paste("TSNE.pdf",sep="_"))
DimPlot(object = Combined, reduction = "tsne", group.by = "condition",cols = color_cond)
DimPlot(object = Combined, reduction = "tsne", group.by = "Cluster", label = TRUE,cols = color_clust)
DimPlot(object = Combined, reduction = "tsne", group.by = "Cell_Type", label = TRUE,cols = color_cells)
dev.off()


setwd("../")

saveRDS(Combined,"3.Combined.rds")

