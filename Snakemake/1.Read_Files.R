#639ee45f
source("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/Snakemake/1.Functions.R")
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
NewDir <- paste0("Rescue_MBSYN_",tool,"_elbow_",elbow,"_remove_mt_",remove_mt,"_remove_ribsomal_",remove_ribsomal)
dir.create(NewDir)
setwd(NewDir)



# ==============================================================================================
# ================================ Setup the Seurat objects ====================================
# ==============================================================================================
# ======== Perform an integrated analysis ====
dir.create("1.QC")
setwd("1.QC")
#debugonce(create_cds)
Return_fun <- create_cds(list_of_files=list_of_files,
                         condition_names=condition_names,
                         min.features =200,min.cells=5,remove_mt=remove_mt,data_10x=data_10x, elbow = elbow,tool=tool,imputation=imputation)
Combined  <- Return_fun$Combined
Data_List <- Return_fun$Data_List
setwd("../")

Mit.index <- grep(pattern = "^MT-|^MT\\.", x = rownames(object), value = FALSE)
rownames(object)[Mit.index]
Mit.index <- grep(pattern = "^MT-|^MT\\.", x = rownames(Combined@assays$RNA@counts), value = FALSE)
Mit.index 


saveRDS(Combined,"1.Combined.rds")
