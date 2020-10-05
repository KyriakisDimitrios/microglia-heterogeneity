
source("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/Snakemake/4.Functions.R")

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





Combined <- readRDS("3.Combined.rds")










# =================================== DF Conditions =================================================
dir.create("4.DF_Conditions")
setwd("4.DF_Conditions")
return_fun_cond <- df_genes(Combined,"condition",top_n=15,logfc.threshold=0, n_cores=4,
    latent.vars=c("nCount_RNA"),
    only.pos=TRUE)
setwd("../")
# ------------------------------------------------------------------------------------------------




# =================================== DF Clusters =================================================
dir.create("4.DF_Clusters")
setwd("4.DF_Clusters")
return_fun_cl <- df_genes(Combined,"Cluster",top_n=15,
    logfc.threshold=0.0, n_cores=4,
    latent.vars=c("nCount_RNA"),only.pos=TRUE)


annotated_heat(Combined,c(1),return_fun_cl$top_markers,"cluster","cluster",ordering="Cluster",Colv=NA,Rowv=NA,One_annot=FALSE)


library(stackoverflow)
cl_tops <- chunk2(return_fun_cl$top_markers,length(unique(Combined$Cluster)))
list_bar_cl <- c()
for (i in cl_tops){
	list_bar_cl <- c(list_bar_cl,i[1:3])
}


gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","MEG3","CALD1")
plot_bar_cells(list_bar_cl[1:6],object=Combined,save=TRUE,title="1-6")
plot_bar_cells(list_bar_cl[7:12],object=Combined,save=TRUE,title="7-12")
plot_bar_cells(list_bar_cl[13:18],object=Combined,save=TRUE,title="13-18")
plot_bar_cells(list_bar_cl[19:24],object=Combined,save=TRUE,title="19-24")
plot_bar_cells(list_bar_cl[25:27],object=Combined,save=TRUE,title="25-27")



pdf('DF_VLN_1.pdf',height=24)
DefaultAssay(Combined)<-"RNA"
plots <-VlnPlot(object = Combined,features = list_bar_cl[1:9],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Combined)<-"integrated"
dev.off()

pdf('DF_VLN_2.pdf',height=24)
DefaultAssay(Combined)<-"RNA"
plots <-VlnPlot(object = Combined,features = list_bar_cl[10:19],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Combined)<-"integrated"
dev.off()

pdf('DF_VLN_3.pdf',height=24)
DefaultAssay(Combined)<-"RNA"
plots <-VlnPlot(object = Combined,features = list_bar_cl[20:27],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Combined)<-"integrated"
dev.off()




setwd("../")