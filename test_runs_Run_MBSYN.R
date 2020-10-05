lapply(list.files("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/icswrapper/R",full.names  = T),source)

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
library(monocle)
# library(scran)
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



# color_cells <-primary.colors(15, steps = 3, no.white = TRUE)



# ================================ SETTING UP ======================================== #
# Number of cells to use
imputation = FALSE
remove_mt=TRUE
remove_ribsomal=TRUE
n_cores=5
elbow = FALSE


# setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-11-20_seurat_elbow_FALSE")
# Combined <- readRDS("Cell_Types_2019-11-20_seurat_elbow_FALSE.rds")
setwd("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/")

NewDir <- paste0(Sys.Date(),"_TESTTTESETESTSTES",tool,"_elbow_",elbow,"_remove_ribsomal_",remove_ribsomal)

dir.create(NewDir)
setwd(NewDir)


# setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-11-13_seurat_elbow_FALSE_SAVER")

# ==============================================================================================
# ================================ Setup the Seurat objects ====================================
# ==============================================================================================
# ======== Perform an integrated analysis ====
dir.create("QC")
setwd("QC")
# debugonce(load_files)
Return_fun <- create_cds2(list_of_files=list_of_files,
                         condition_names=condition_names,
                         min.features =200,min.cells=5,
                         remove_mt=TRUE,
                         data_10x=data_10x,
                          elbow = elbow,
                         tool=tool,
                         imputation=imputation,
                        n_cores=2)

Combined  <- Return_fun$Combined
Data_List <- Return_fun$Data_List
setwd("../")

# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------

#saveRDS(Combined,paste0("Merged_",NewDir,".rds"))
# setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-11-13_seurat_elbow_FALSE_SAVER")
# Combined <- readRDS("Clustered_2019-12-04_seurat_elbow_FALSE_remove_ribsomal_TRUEVALIDATION.rds")

# ==============================================================================================
# ======================================= Clustering  ==========================================
# ==============================================================================================
dir.create("Clustering")
setwd("Clustering")
# Combined <- ReduceDim(Combined,method="umap",project=project)$Combined
# debugonce(reduce_dim)
Combined <- reduce_dim(Combined,project=project)$Combined

pdf(paste(Sys.Date(),project,"tsne","Conditions.pdf",sep="_"))
plot_cells(Combined,target="condition",leg_pos="right",save=TRUE,ncol=1)
dev.off()
pdf(paste(Sys.Date(),project,"tsne","Cluster.pdf",sep="_"))
# debugonce(plot_cells)
plot_cells(Combined,target="Cluster",leg_pos="right",save=TRUE,ncol=1)
dev.off()


plot_nFeatures(Combined,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)
plot_tot_mRNA(Combined,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)


if(tolower(tool)=="seurat" & elbow){
    p3 <- DimPlot(object = Combined, reduction = "umap", group.by = "condition",cols = color_cond)
    p4 <- DimPlot(object = Combined, reduction = "umap", group.by = "Cluster", label = TRUE,cols = color_clust)

    pdf(paste(Sys.Date(),project,"umap","Seurat.pdf",sep="_"))
    print(p3)
    print(p4)
    dev.off()
    p3 <- DimPlot(object = Combined, reduction = "tsne", group.by = "condition",cols = color_cond)
    p4 <- DimPlot(object = Combined, reduction = "tsne", group.by = "Cluster", label = TRUE,cols = color_clust)

    pdf(paste(Sys.Date(),project,"tsne","Seurat.pdf",sep="_"))
    print(p3)
    print(p4)
    dev.off()
}

# ------------------------------------------------------------------------------------------------

# ====================== BARPLOT CELL BELONGS TO WHICH CLUSTERS ====================================
data <- data.frame(cbind(factor(Combined$Cluster),paste(" ",Combined$condition)))
colnames(data)<- c("Cluster","Condition")
pdf(paste(Sys.Date(),'Barplot_num-Cond_per_Cluster.pdf',sep="_"))
print(ggplot(data, aes(Cluster)) +  geom_bar(aes(fill = Condition))+
          guides(col=guide_legend(ncol=2,))+
          scale_fill_manual(values =color_list$condition)+
          theme(legend.position="bottom") )
dev.off()
setwd("../")
# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------

saveRDS(Combined,paste0("Clustered_",NewDir,".rds"))

# =================================== DF Conditions =================================================
dir.create("DF_Conditions")
setwd("DF_Conditions")
# debugonce(df_genes)
return_fun_cond <- df_genes(Combined,"condition",top_n=15,logfc.threshold=0, n_cores=4,
    latent.vars=c("nCount_RNA"),
    only.pos=TRUE,color_list = color_list)
setwd("../")



# ------------------------------------------------------------------------------------------------




# =================================== DF Clusters =================================================
dir.create("DF_Clusters")
setwd("DF_Clusters")
return_fun_cl <- df_genes(Combined,"Cluster",top_n=15,
    logfc.threshold=0.0, n_cores=4,
    latent.vars=c("nCount_RNA"),only.pos=TRUE,color_list = color_list)


annotated_heat(Combined,c(1),return_fun_cl$top_markers,"cluster","cluster",ordering="Cluster",Colv=NA,Rowv=NA,One_annot=FALSE,color_list=NULL)


library(stackoverflow)
cl_tops <- chunk2(return_fun_cl$top_markers,length(unique(Combined$Cluster)))
list_bar_cl <- c()
for (i in cl_tops){
	list_bar_cl <- c(list_bar_cl,i[1:3])
}


library(clusterProfiler)
dir.create("GSEA")
setwd("GSEA")
# debugonce(gene_set_enrich)
sub_degs <- return_fun_cl$degs[return_fun_cl$degs$p_val_adj<0.05 & return_fun_cl$degs$avg_logFC>1.5,]
universe <- rownames(Combined@assays$RNA@counts)

symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
test <- split(symbols,sub_degs$cluster)

x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,ont="BP")
pdf("enrichGO.pdf",width=12,height=10)
dotplot(x, showCategory=5, includeAll=FALSE)
dev.off()

x2=compareCluster(test, fun='groupGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",ont="BP",level=2)

pdf("groupGO.pdf",width=12)
dotplot(x2, showCategory=10, includeAll=FALSE)
dev.off()

list_symbols <- list()
iter=0
for (i in names(test)){
	iter=iter+1
	eg=bitr(as.vector(test[[i]]),fromType = "SYMBOL",toType=c("SYMBOL","ENTREZID"),OrgDb='org.Mm.eg.db')$ENTREZID
	list_symbols[[i]] <-c(eg)
}

x3=compareCluster(list_symbols, fun='enrichKEGG', organism="mmu",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05)

pdf("enrichKEGG.pdf",width=12)
dotplot(x3, showCategory=10, includeAll=FALSE)
dev.off()




for (cl in unique(sub_degs$cluster)){
	dir.create(paste0("Cluster_",cl))
	setwd(paste0("Cluster_",cl))
	gene_symbols <- sub_degs[sub_degs$cluster==cl,]$gene
	ego <-gene_set_enrich(gene_symbols,organism=organism,ontology=c("BP","MF"),qval_thres=0.05,save=TRUE,title=cl)
	# pdf(paste(Sys.Date(),"GSEA_Dotplot","Cluster",cl,".pdf",sep="_"),height=14,width=18)
	# print(barplot(ego, showCategory=30))
	# dev.off()
	# debugonce(david_enrich)
	das <- david_enrich(gene_symbols=gene_symbols,
		idents=NULL,organism=organism,
		qval_thres=0.05,save=TRUE,title=paste0("CL",cl))
	setwd("../")
}


setwd("../")
setwd("../")

# ------------------------------------------------------------------------------------------------






dir.create("Cell_Types_Garnett")
setwd("Cell_Types_Garnett")
library(garnett)
Combined <- garnett_assignment(Combined ,organism="mouse",check_markers=FALSE)
# Combined <- garnett_assignment(Combined ,organism="skata",check_markers=FALSE)
pdf("Garnett_Cell_type.pdf")
plot_cells(Combined,target="cell_type")
dev.off()
pdf("Garnett_cluster_ext_type.pdf")
plot_cells(Combined,target="cluster_ext_type")
dev.off()
pdf("Garnett_garnett_cluster.pdf")
plot_cells(Combined,target="garnett_cluster")
dev.off()
setwd("../")




# =================================== CELL TYPES =================================================
dir.create("Cell_Types")
setwd("Cell_Types")
gene_list_name <- "Cell_Type"
file <- paste0(WORKDIR,"/Gene_Lists/Paper_Cell_types.txt")
color_list$Cell_Type <- color_list$Cell_Type[-4]

res_cell_type_assignment <- cell_type_assignment(object=Combined,tab_name="Cell_Type",file=file)#,color_list = color_list)

Combined = res_cell_type_assignment$object
r_annot = res_cell_type_assignment$r_annot

# ============ ADVANCED QUALITY FILTERING =========================
cell_types <- Combined$Cell_Type
cell_types <- as.vector(Combined$Cluster)

save_cell_types<-Combined$Cell_Type
cell_types[ which(Combined$Cluster %in% c(1))] <- "Astrocytes"
cell_types[ which(Combined$Cluster %in% c(2))] <- "Microglia"
cell_types[ which(Combined$Cluster %in% c(3))] <- "Oligodendrocytes"
cell_types[ which(Combined$Cluster %in% c(4))] <- "Endothelial"
cell_types[ which(Combined$Cluster %in% c(6,7))] <- "Ependymal"
cell_types[ which(Combined$Cluster %in% c(5,8))] <- "Neurons"
cell_types[ which(Combined$Cluster %in% c(9))] <- "Hybrid"
Combined$Cell_Type <-cell_types

# ================================================================


gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","MEG3","CALD1","MYL9","VTN")
plot_bar_cells(gene_list,object=Combined,save=TRUE)
# debugonce(pie_per_CellType)
# library(d)
pie_per_CellType(Combined,target="Cell_Type")

setwd("../")



# ======================================== DF CELL TYPES  =================================
# ======================================== DF CELL TYPES  =================================
# ======================================== DF CELL TYPES  =================================
dir.create("DF_Cell_Types")
setwd("DF_Cell_Types")
return_fun_cell_types <- df_genes(Combined,"Cell_Type",top_n=50,logfc.threshold=1, n_cores=4,latent.vars=c("nCount_RNA"))
gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","MEG3","CALD1")
plot_bar_cells(gene_list,object=Combined,save=TRUE)
pie_per_CellType(Combined,target="Cell_Type")
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
# ------------------------------------------------------------------------------------------------

plot_bar_cells(list_bar_cl[1:6],object=Combined,save=TRUE,title="1-6",target="Cluster")
plot_bar_cells(list_bar_cl[7:12],object=Combined,save=TRUE,title="7-12",target="Cluster")
plot_bar_cells(list_bar_cl[13:18],object=Combined,save=TRUE,title="13-18",target="Cluster")
plot_bar_cells(list_bar_cl[19:24],object=Combined,save=TRUE,title="19-24",target="Cluster")
plot_bar_cells(list_bar_cl[25:27],object=Combined,save=TRUE,title="25-27",target="Cluster")

saveRDS(Combined,paste0("Cell_Types_",NewDir,".rds"))










cell_types <- Combined$Cell_Type
cell_types <- as.vector(Combined$Cluster)
save_cell_types<-Combined$Cell_Type
cell_types[ which(Combined$Cluster %in% c(1))] <- "Astrocytes"
cell_types[ which(Combined$Cluster %in% c(2))] <- "Microglia"
cell_types[ which(Combined$Cluster %in% c(3))] <- "Oligodendrocytes"
cell_types[ which(Combined$Cluster %in% c(4))] <- "Endothelial"
cell_types[ which(Combined$Cluster %in% c(5))] <- "Hybrid 1"
cell_types[ which(Combined$Cluster %in% c(6))] <- "Ependymal"
cell_types[ which(Combined$Cluster %in% c(7))] <- "Epithelial"
cell_types[ which(Combined$Cluster %in% c(8))] <- "Neurons"
cell_types[ which(Combined$Cluster %in% c(9))] <- "Hybrid 2"
Combined$Cell_Type <-cell_types

color_clust_ord_cl <- color_clust[c(1,4,6,7,5,9,2,8,3)]
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_clust_ord_cl,State=color_clust)




p1 <- plot_cells(Combined,target="condition",leg_pos="right",save=FALSE,ncol=1,color_list = color_list,pt.size=2)
p2 <-plot_cells(Combined,target="Cluster",leg_pos="right",save=FALSE,ncol=1,color_list = color_list,pt.size=2)
p3 <-plot_cells(Combined,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,color_list = color_list,pt.size=2)

pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

pdf(paste(Sys.Date(),project,"Cell_Type_tsne.pdf",sep="_"))
print(p3)
dev.off()
pdf(paste(Sys.Date(),project,"Condition_tsne.pdf",sep="_"))
print(p1)
dev.off()
pdf(paste(Sys.Date(),project,"Cluster_tsne.pdf",sep="_"))
print(p2)
dev.off()


p1 <- plot_cells(Combined,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list = color_list,pt.size=2)
p2 <-plot_cells(Combined,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list = color_list,pt.size=2)
p3 <-plot_cells(Combined,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list = color_list,pt.size=2)

pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

pdf(paste(Sys.Date(),project,"Cell_Type_umap.pdf",sep="_"))
print(p3)
dev.off()
pdf(paste(Sys.Date(),project,"Condition_umap.pdf",sep="_"))
print(p1)
dev.off()
pdf(paste(Sys.Date(),project,"Cluster_umap.pdf",sep="_"))
print(p2)
dev.off()
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------





# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ///////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ///////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ///////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
library(ICSWrapper)
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




# ================================ SETTING UP ======================================== #
# Number of cells to use
imputation = FALSE
remove_mt=TRUE
remove_ribsomal=TRUE
n_cores=5
elbow = FALSE
setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-12-11_seurat_elbow_FALSE_remove_ribsomal_TRUE")
# Combined <- readRDS("Cell_Types_2019-12-05_seurat_elbow_FALSE_remove_ribsomal_TRUE.rds")
load("Paper_WORKSPACE.RData")
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ///////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ///////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ///////////////// <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Combined$Cell_Type <- Combined$cluster_ext_type


save.image("Paper_WORKSPACE.RData")









# =============================== Microglia Analysis Integraded  ===============================
# =============================== Microglia Analysis Integraded  ===============================
# ======================================== Microglia Analysis  =================================
dir.create("Microglia_Analysis")
setwd("Microglia_Analysis")
color_clust_microglia <- c("#42858C" ,rgb(128, 128, 192, maxColorValue = 255),
    rgb(255,128,192, maxColorValue = 255),
    rgb(128,0,64, maxColorValue = 255),
    rgb(128,0,255, maxColorValue = 255),"blue")
color_cells_microglia <- c(rgb(128,0,128, maxColorValue = 255))
color_list_microglia <- list(condition=color_cond,Cluster=color_clust_microglia,Cell_Type=color_cells_microglia,State=color_clust_microglia)




debugonce(object_subset)
Micro_Sub <- object_subset(object=Combined,Conditions=NULL,Clusters=NULL,Cell_Types=c("Microglia") )
DefaultAssay(Micro_Sub) <- "integrated"
Micro_Sub <- reduce_dim(Micro_Sub,project=project,resolution=c(1.2))$Combined



gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster")

DimPlot(Micro_Sub,group.by = "condition",reduction = "umap",cols=color_list_microglia$condition)
DimPlot(Micro_Sub,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)
DimPlot(Micro_Sub,group.by = "condition",reduction = "tsne",cols=color_list_microglia$condition)
DimPlot(Micro_Sub,group.by = "Cluster",reduction = "tsne",cols=color_list_microglia$Cluster)

p1<-plot_cells(Micro_Sub,target="condition",leg_pos="bottom",save=FALSE,ncol=1,reduction="umap",color_list=color_list_microglia,pt.size=2)
p2<-plot_cells(Micro_Sub,target="Cluster",leg_pos="bottom",save=FALSE,ncol=1,reduction="umap",color_list=color_list_microglia,pt.size=2)
p3<-plot_cells(Micro_Sub,target="Cell_Type",leg_pos="bottom",save=FALSE,ncol=1,reduction="umap",color_list=color_list_microglia,pt.size=2)

pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=8,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
print(p1+theme(legend.position="bottom"))
print(p2+theme(legend.position="bottom"))
print(p3+theme(legend.position="bottom"))
dev.off()

p1<-plot_cells(Micro_Sub,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list=color_list_microglia,pt.size=2)
p2<-plot_cells(Micro_Sub,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list=color_list_microglia,pt.size=2)
p3<-plot_cells(Micro_Sub,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list=color_list_microglia,pt.size=2)

pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()


pdf("Microglia2.pdf")
DimPlot(object = Micro_Sub2, reduction = "umap", group.by = "condition",cols = color_list_microglia[["condition"]],pt.size=2)+theme(legend.position = "bottom")
DimPlot(object = Micro_Sub2, reduction = "umap", group.by = "Cluster",cols = color_list_microglia[["Cluster"]],pt.size=2)+theme(legend.position = "bottom")
dev.off()


# =================================== DF Conditions =================================================
# =================================== DF Conditions =================================================
dir.create("DF_Conditions")
setwd("DF_Conditions")
return_fun_cond_micro <- df_genes(Micro_Sub,"condition",top_n=15,logfc.threshold=0,
    n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE,color_list=color_list_microglia)
setwd("../")
# ---------------------------------------------------------------------------------------------




# =================================== DF Clusters =================================================
# =================================== DF Clusters =================================================
dir.create("DF_Clusters")
setwd("DF_Clusters")
return_fun_cl_micro <- df_genes(Micro_Sub,"Cluster",top_n=30,logfc.threshold=0 ,
    n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE,color_list=color_list_microglia)

gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster",title="Homeostatic")


library(clusterProfiler)
dir.create("GSEA")
setwd("GSEA")
# debugonce(gene_set_enrich)
sub_degs <- return_fun_cl_micro$degs[return_fun_cl_micro$degs$p_val_adj<0.05 & return_fun_cl_micro$degs$avg_logFC>1.5,]
universe <- rownames(Micro_Sub@assays$RNA@counts)

symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
test <- split(symbols,sub_degs$cluster)

x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,ont="BP")
pdf("enrichGO.pdf",width=12,height=10)
dotplot(x, showCategory=5, includeAll=FALSE)
dev.off()
setwd("../")
setwd("../")









# ============================= REANALUSIS ===============================
# ============================= REANALUSIS ===============================
# ============================= REANALUSIS ===============================

dir.create("Reanalysis")
setwd("Reanalysis")
Micro_Sub2 <- object_subset(object=Micro_Sub,Conditions=NULL,Clusters=c(1,2,3),Cell_Types=NULL)
DefaultAssay(Micro_Sub2) <- "RNA"
Micro_Sub2 <- reduce_dim(Micro_Sub2,project=project,resolution=c(1.1))$Combined

p1<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list=color_list_microglia,pt.size=2)
p2<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list=color_list_microglia,pt.size=2)
p3<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="umap",color_list=color_list_microglia,pt.size=2)

pdf(paste(Sys.Date(),project,"_Microglia2_umap.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

p1<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list=color_list_microglia,pt.size=2)
p2<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list=color_list_microglia,pt.size=2)
p3<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",color_list=color_list_microglia,pt.size=2)

pdf(paste(Sys.Date(),project,"_Microglia2_tsne.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

# =================================== DF Conditions =================================================
# =================================== DF Conditions =================================================
dir.create("DF_Conditions")
setwd("DF_Conditions")
return_fun_cond_micro2 <- df_genes(Micro_Sub2,"condition",top_n=15,logfc.threshold=0,
    n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE,color_list=color_list_microglia)

micro2_all_cond_df <- return_fun_cond_micro2$degs$gene[return_fun_cond_micro2$degs$p_val<0.05]
annotated_heat(Micro_Sub2,c(1),micro2_all_cond_df,"condition","condition",
    ordering="condition",Colv=NA,Rowv=NA,One_annot=FALSE,color_list=color_list_microglia)

setwd("../")
# ---------------------------------------------------------------------------------------------


# =================================== DF Clusters =================================================
# =================================== DF Clusters =================================================
dir.create("DF_Clusters")
setwd("DF_Clusters")
return_fun_cl_micro2 <- df_genes(Micro_Sub2,"Cluster",top_n=30,
    logfc.threshold=0 ,n_cores=4,latent.vars=c("nCount_RNA"),
    only.pos=TRUE,color_list=color_list_microglia)
color_list <- color_list_microglia
gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster",title="Homeostatic")

library(clusterProfiler)
dir.create("GSEA")
setwd("GSEA")
# debugonce(gene_set_enrich)
sub_degs <- return_fun_cl_micro2$degs[return_fun_cl_micro2$degs$p_val_adj<0.05 & return_fun_cl_micro2$degs$avg_logFC>0.5,]
universe <- rownames(Micro_Sub2@assays$RNA@counts)

symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
test <- split(symbols,sub_degs$cluster)

x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,ont="BP")
pdf("enrichGO.pdf",width=12,height=10)
dotplot(x, showCategory=5, includeAll=FALSE)
dev.off()
setwd("../")

pdf('DF_VLN_Seurat.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = return_fun_cl_micro2$top_markers,group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()


dir.create("monocle")
setwd("monocle")
micro_mat  <- as.matrix(Micro_Sub2@assays$RNA@counts)
micro_cl <- Micro_Sub2$Cluster
object <- seurat_to_monocle(Micro_Sub2)
target="Cluster"
all_diff_test_res_condition <- differentialGeneTest(object, fullModelFormulaStr=paste("~",target,sep=""),cores=n_cores, reducedModelFormulaStr ="~nCount_RNA")
# == Reduce Matrix (remove not significant genes)
all_diff_test_res_condition <- all_diff_test_res_condition %>% arrange(qval)
monocle_helper <- df_genes_logfc(object,target,signif_genes=all_diff_test_res_condition,top_n=50,qval_thres=0.05,fc_thres=0,each_cl=TRUE,n_cores=4)
pbmc.markers2 <- monocle_helper$df_pval_genes


top <- pbmc.markers2[pbmc.markers2$qval<0.05,]

# top10<-top%>%arrange(qval)
# top10_genes <-top10$gene_short_name[1:50]

top10 <- top %>% group_by(cluster) %>% top_n(n = 15, wt = -qval)
top10_genes<- unique(top10$gene_short_name)

annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_Microglia",Rowv=FALSE,
                  title="DF_genes_Pvalue",
                  ordering="Cluster")
annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_Microglia",Rowv=TRUE,
                  title="DF_genes_CL",
                  ordering="Cluster")

top10<-top%>%arrange(-avg_logFC)
top10_genes <-top10$gene_short_name[1:50]

top10 <- top %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
top10_genes<- unique(top10$gene_short_name)
annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_Microglia",Rowv=FALSE,
                  title="DF_genes_FC",
                  ordering="Cluster")

annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top$gene_short_name,
                  gene_list_name="DF_Microglia",Rowv=FALSE,
                  title="DF_genes_FC",
                  ordering="Cluster")


write.table(pbmc.markers2,"Monocle_Top_qvalue.tsv",sep="\t")

pdf('DF_VLN_Monocle.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = top10_genes,group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
plots
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()


setwd("../") # End Monocle























# # ======================================== Microglia Analysis  =================================
# dir.create("Microglia_Analysis_counts")
# setwd("Microglia_Analysis_counts")
# # debugonce(object_subset)
# Micro_Sub <- object_subset(object=Combined,Conditions=NULL,Clusters=NULL,Cell_Types=c("Microglia") )
# DefaultAssay(Micro_Sub) <- "RNA"
# Micro_Sub <- reduce_dim(Micro_Sub,project=project,resolution=c(1.2))$Combined



# gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster")




# p1<-plot_cells(Micro_Sub,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="umap")
# p2<-plot_cells(Micro_Sub,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="umap")
# p3<-plot_cells(Micro_Sub,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="umap")

# pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# p1<-plot_cells(Micro_Sub,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="tsne")
# p2<-plot_cells(Micro_Sub,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="tsne")
# p3<-plot_cells(Micro_Sub,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="tsne")

# pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# # =================================== DF Conditions =================================================
# # =================================== DF Conditions =================================================
# dir.create("DF_Conditions")
# setwd("DF_Conditions")
# return_fun <- df_genes(Micro_Sub,"condition",top_n=15,logfc.threshold=0.15,n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE)
# setwd("../")
# # ---------------------------------------------------------------------------------------------




# # =================================== DF Clusters =================================================
# # =================================== DF Clusters =================================================
# dir.create("DF_Clusters")
# setwd("DF_Clusters")
# return_fun <- df_genes(Micro_Sub,"Cluster",top_n=30,logfc.threshold=0 ,n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE)


# gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster",title="Homeostatic")


# library(clusterProfiler)
# dir.create("GSEA")
# setwd("GSEA")
# # debugonce(gene_set_enrich)
# sub_degs <- return_fun$degs[return_fun$degs$p_val_adj<0.05 & return_fun$degs$avg_logFC>1.5,]
# universe <- rownames(Micro_Sub@assays$RNA@counts)

# symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
# test <- split(symbols,sub_degs$cluster)

# x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05,ont="BP")
# pdf("enrichGO.pdf",width=12,height=10)
# dotplot(x, showCategory=5, includeAll=FALSE)
# dev.off()
# setwd("../")
# setwd("../")

# # =============================== PAIRWISE DF ===============================================
# # =============================== PAIRWISE DF ===============================================
# # dir.create("DF_Pairwise")
# # setwd("DF_Pairwise")
# # Idents(Micro_Sub) <- Micro_Sub$Cluster
# # cl_combinations <- combn(levels(Micro_Sub$Cluster),2)

# pairwise_df <- function (comb,cl_combinations,object){
#     title <- paste(cl_combinations[,comb],collapse = "_")
#     dir.create(title)
#     setwd(title)
#     target <- "Cluster"
#     idents <- as.vector(cl_combinations[,comb])
#     ident.1 <- idents[1]
#     ident.2 <- idents[2]


#     cells_index <- object$Cluster%in%c(ident.1,ident.2)
#     temp_object <- subset(x=object,cells=colnames(object)[cells_index])


#     temp_object <- ScaleData(temp_object)
#     DefaultAssay(object) <- "RNA"
#     pbmc.markers <- FindAllMarkers(object = temp_object,
#                                            assay ="RNA",
#                                            logfc.threshold=0,
#                                            only.pos = TRUE,
#                                            test.use = "MAST",latent.vars = c("nCount_RNA"))
#     pbmc.markers$gene <- rownames(pbmc.markers)
#     qvalue <- p.adjust(pbmc.markers$p_val, method = "BH",n=dim(temp_object@assays$RNA@counts)[1])
#     pbmc.markers$qvalue <- qvalue
#     top <- pbmc.markers[pbmc.markers$qvalue<0.05,]
#     top10 <- top %>% top_n(n = 50, wt = abs(avg_logFC))
#     top10_genes<- top10$gene
#     temp <- temp_object[,temp_object$Cluster%in%c(ident.1,ident.2)]
#     temp$Cluster <- as.factor(as.vector(temp$Cluster))

#     # debugonce(annotated_heat)
#     annotated_heat(object=temp,
#                       row_annotation=c(1),
#                       gene_list=top10_genes,
#                       Rowv=TRUE,
#                       gene_list_name="DF_genes",
#                       title=title,
#                       ordering="Cluster")

#     DefaultAssay(temp_object) <- "integrated"
#     write.table(pbmc.markers, file = paste0(Sys.Date(),"_TO_EXP_each_",target,"_",title,".tsv"),row.names=FALSE, na="", sep="\t")

#     test <- split(top$gene,top$cluster)

#     x=compareCluster(test, fun='enrichGO', OrgDb='org.Hs.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                             pvalueCutoff  = 0.05,ont="BP")
#     pdf(paste0(Sys.Date(),"_enrichGO_BP_",target,"_",title,".pdf"),width=12,height=10)
#     dotplot(x, showCategory=5, includeAll=FALSE)
#     dev.off()

#     x=compareCluster(test, fun='enrichGO', OrgDb='org.Hs.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                             pvalueCutoff  = 0.05,ont="MF")
#     pdf(paste0(Sys.Date(),"_enrichGO_MF_",target,"_",title,".pdf"),width=12,height=10)
#     dotplot(x, showCategory=5, includeAll=FALSE)
#     dev.off()
#     pbmc.markers[c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3"),]
#     setwd("../")
# }


# # lapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,cl_combinations=cl_combinations,object=Micro_Sub)

# # Homeostatic <- c(grep("FCRLS",return_fun$degs$gene),
# # grep("LGMN",return_fun$degs$gene),
# # grep("MARCKS",return_fun$degs$gene),
# # grep("CSF1R",return_fun$degs$gene),
# # grep("P2RY12",return_fun$degs$gene),
# # grep("SERINC3",return_fun$degs$gene))

# # return_fun$degs[Homeostatic,]
# # setwd("../")


# # ============================= REANALUSIS ===============================
# # ============================= REANALUSIS ===============================
# # ============================= REANALUSIS ===============================

# dir.create("Reanalysis")
# setwd("Reanalysis")
# Micro_Sub2 <- object_subset(object=Micro_Sub,Conditions=NULL,Clusters=c(1,2,3,4,6),Cell_Types=NULL )
# DefaultAssay(Micro_Sub2) <- "RNA"
# Micro_Sub2 <- reduce_dim(Micro_Sub2,project=project,resolution=c(1.2))$Combined
# gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
# plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster")




# p1<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="umap",pt.size=2)
# p2<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="umap",pt.size=2)
# p3<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="umap",pt.size=2)

# pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# p1<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",pt.size=2)
# p2<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",pt.size=2)
# p3<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,reduction="tsne",pt.size=2)

# pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# # =================================== DF Conditions =================================================
# # =================================== DF Conditions =================================================
# dir.create("DF_Conditions")
# setwd("DF_Conditions")
# return_fun <- df_genes(Micro_Sub2,"condition",top_n=15,logfc.threshold=0.15,n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE)
# setwd("../")
# # ---------------------------------------------------------------------------------------------

# # =================================== DF Clusters =================================================
# # =================================== DF Clusters =================================================
# dir.create("DF_Clusters")
# setwd("DF_Clusters")
# return_fun <- df_genes(Micro_Sub2,"Cluster",top_n=30,logfc.threshold=0 ,n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE)


# gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
# plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster",title="Homeostatic")


# library(clusterProfiler)
# dir.create("GSEA")
# setwd("GSEA")
# # debugonce(gene_set_enrich)
# sub_degs <- return_fun$degs[return_fun$degs$p_val_adj<0.05 & return_fun$degs$avg_logFC>1.5,]
# universe <- rownames(Micro_Sub2@assays$RNA@counts)

# symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
# test <- split(symbols,sub_degs$cluster)

# x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05,ont="BP")
# pdf("enrichGO.pdf",width=12,height=10)
# dotplot(x, showCategory=5, includeAll=FALSE)
# dev.off()
# setwd("../")
# setwd("../")
# setwd("../")

# =============================== PAIRWISE DF ===============================================
# =============================== PAIRWISE DF ===============================================
# dir.create("DF_Pairwise")
# setwd("DF_Pairwise")
# Idents(Micro_Sub2) <- Micro_Sub2$Cluster
# cl_combinations <- combn(levels(Micro_Sub2$Cluster),2)

# lapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,cl_combinations=cl_combinations,object=Micro_Sub2)

# Homeostatic <- c(grep("FCRLS",return_fun$degs$gene),
# grep("LGMN",return_fun$degs$gene),
# grep("MARCKS",return_fun$degs$gene),
# grep("CSF1R",return_fun$degs$gene),
# grep("P2RY12",return_fun$degs$gene),
# grep("SERINC3",return_fun$degs$gene))

# return_fun$degs[Homeostatic,]
# setwd("../")

# setwd("../")
# ----------------------------------------------------------------------------------------------



#
# =============================== PAIRWISE DF ===============================================
# =============================== PAIRWISE DF ===============================================
# dir.create("DF_Pairwise")
# setwd("DF_Pairwise")
# Idents(Micro_Sub) <- Micro_Sub$Cluster
# cl_combinations <- combn(levels(Micro_Sub$Cluster),2)


# lapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,cl_combinations=cl_combinations,object=Micro_Sub)

# Homeostatic <- c(grep("FCRLS",return_fun$degs$gene),
# grep("LGMN",return_fun$degs$gene),
# grep("MARCKS",return_fun$degs$gene),
# grep("CSF1R",return_fun$degs$gene),
# grep("P2RY12",return_fun$degs$gene),
# grep("SERINC3",return_fun$degs$gene))

# return_fun$degs[Homeostatic,]
# setwd("../")

pdf('VLN_APOE.pdf')
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = c("APOE"),group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
dev.off()

setwd("../")


dir.create("Reanalysis_Counts")
setwd("Reanalysis_Counts")
Micro_Sub2 <- object_subset(object=Micro_Sub,Conditions=NULL,Clusters=c(1,2,3),Cell_Types=NULL)
DefaultAssay(Micro_Sub2) <- "RNA"
Micro_Sub2 <- reduce_dim(Micro_Sub2,project=project,resolution=c(1.2))$Combined

p1<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list_microglia)
p2<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list_microglia)
p3<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list_microglia)

pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

p1<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=TRUE,ncol=1,reduction="tsne",color_list=color_list_microglia)
p2<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=TRUE,ncol=1,reduction="tsne",color_list=color_list_microglia)
p3<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,reduction="tsne",color_list=color_list_microglia)

pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

# =================================== DF Conditions =================================================
# =================================== DF Conditions =================================================
dir.create("DF_Conditions")
setwd("DF_Conditions")
return_fun_cond_micro2 <- df_genes(Micro_Sub2,"condition",top_n=15,logfc.threshold=0,
    n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE,color_list=color_list_microglia)


setwd("../")
# ---------------------------------------------------------------------------------------------

# =================================== DF Clusters =================================================
# =================================== DF Clusters =================================================
dir.create("DF_Clusters")
setwd("DF_Clusters")
return_fun_cl_micro2 <- df_genes(Micro_Sub2,"Cluster",top_n=30,
    logfc.threshold=0 ,n_cores=4,latent.vars=c("nCount_RNA"),
    only.pos=TRUE,color_list=color_list_microglia)
color_list_g <- color_list
color_list <- color_list_microglia
color_list$Cell_Type <- color_list$Cluster
gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster",title="Homeostatic")

library(clusterProfiler)
dir.create("GSEA")
setwd("GSEA")
sub_degs <- return_fun_cl_micro2$degs[return_fun_cl_micro2$degs$p_val_adj<0.05 & return_fun_cl_micro2$degs$avg_logFC>0.5,]
universe <- rownames(Micro_Sub2@assays$RNA@counts)

symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
test <- split(symbols,sub_degs$cluster)

x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
                        pvalueCutoff  = 0.05,ont="BP")
pdf("enrichGO.pdf",width=12,height=10)
dotplot(x, showCategory=5, includeAll=FALSE)
dev.off()
setwd("../")


pdf('DF_VLN_Seurat.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = return_fun_cl_micro2$top_markers,group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()

dir.create("monocle")
setwd("monocle")
micro_mat  <- as.matrix(Micro_Sub2@assays$RNA@counts)
micro_cl <- Micro_Sub2$Cluster
object <- seurat_to_monocle(Micro_Sub2)
target="Cluster"
all_diff_test_res_condition <- differentialGeneTest(object, fullModelFormulaStr=paste("~",target,sep=""),cores=n_cores, reducedModelFormulaStr ="~nCount_RNA")
# == Reduce Matrix (remove not significant genes)
all_diff_test_res_condition <- all_diff_test_res_condition %>% arrange(qval)
monocle_helper <- df_genes_logfc(object,target,signif_genes=all_diff_test_res_condition,top_n=50,qval_thres=0.05,fc_thres=0,each_cl=TRUE,n_cores=4)
pbmc.markers2 <- monocle_helper$df_pval_genes


top <- pbmc.markers2[pbmc.markers2$qval<0.05,]

# top10<-top%>%arrange(qval)
# top10_genes <-top10$gene_short_name[1:50]

top10 <- top %>% group_by(cluster) %>% top_n(n = 15, wt = -qval)
top10_genes<- unique(top10$gene_short_name)

annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_Microglia",Rowv=FALSE,
                  title="DF_genes_Pvalue",
                  ordering="Cluster",color_list=color_list_microglia)
annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_Microglia",Rowv=TRUE,
                  title="DF_genes_CL",
                  ordering="Cluster",color_list=color_list_microglia)

top10<-top%>%arrange(-avg_logFC)
top10_genes <-top10$gene_short_name[1:50]

top10 <- top %>% group_by(cluster) %>% top_n(n = 15, wt = avg_logFC)
top10_genes<- unique(top10$gene_short_name)
annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_Microglia",Rowv=FALSE,
                  title="DF_genes_FC",
                  ordering="Cluster")#,color_list=color_list_microglia)


write.table(pbmc.markers2,"Monocle_Top_qvalue.tsv",sep="\t")

pdf('DF_VLN_Monocle.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = top10_genes[1:10],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()


setwd("../") # End Monocle

pdf('VLN_APOE.pdf')
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = c("APOE"),group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
dev.off()


dir.create("Violin_DF")
setwd("Violin_DF")
pdf(paste0('VLN_DF.pdf'))

for (gene in top$gene_short_name){

    DefaultAssay(Micro_Sub2)<-"RNA"
    plots <-VlnPlot(object = Micro_Sub2,features = c(gene),group.by = "Cluster",ncol = 1,combine=FALSE)
    plots <- lapply(
      X = plots,
      FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
    )
    print(plots)
}
dev.off()




setwd("../") # End DF_Clusters




# ======== FILES
dir.create("Files")
setwd("Files")
write.table(rownames(Combined@assays$RNA),"Genes.tsv")
write.table(cbind(colnames(Combined@assays$RNA),as.vector(Combined$condition)),"CellNames_BrainRegion.tsv")
write.table(Combined@reductions$pca@cell.embeddings,"PCA.tsv")
write.table(Combined@reductions$tsne@cell.embeddings,"TSNE.tsv")
write.table(Combined@reductions$umap@cell.embeddings,"UMAP.tsv")
setwd("../")
# =============================== PAIRWISE DF ===============================================
# =============================== PAIRWISE DF ===============================================
# dir.create("DF_Pairwise")
# setwd("DF_Pairwise")
# Idents(Micro_Sub2) <- Micro_Sub2$Cluster
# cl_combinations <- combn(levels(Micro_Sub2$Cluster),2)

# lapply(c(1:dim(cl_combinations)[2]),FUN=pairwise_df,cl_combinations=cl_combinations,object=Micro_Sub2)

# Homeostatic <- c(grep("FCRLS",return_fun$degs$gene),
# grep("LGMN",return_fun$degs$gene),
# grep("MARCKS",return_fun$degs$gene),
# grep("CSF1R",return_fun$degs$gene),
# grep("P2RY12",return_fun$degs$gene),
# grep("SERINC3",return_fun$degs$gene))

# return_fun$degs[Homeostatic,]
# setwd("../")

# setwd("../")




# cl_combinations=cl_combinations,
# object=Micro_Sub2
# title <- paste(cl_combinations[,comb],collapse = "_")
# dir.create(title)
# setwd(title)
# target <- "Cluster"
# idents <- as.vector(cl_combinations[,comb])
# ident.1 <- idents[1]
# ident.2 <- idents[2]


# cells_index <- object$Cluster%in%c(ident.1,ident.2)
# temp_object <- subset(x=object,cells=colnames(object)[cells_index])


# temp_object <- ScaleData(temp_object)
# DefaultAssay(object) <- "RNA"
# pbmc.markers <- FindAllMarkers(object = temp_object,
#                                        assay ="RNA",
#                                        logfc.threshold=0,min.pct = 0,
#                                        only.pos = TRUE,
#                                        test.use = "MAST",latent.vars = c("nCount_RNA"))
# pbmc.markers$gene <- rownames(pbmc.markers)
# qvalue <- p.adjust(pbmc.markers$p_val, method = "BH",n=dim(temp_object@assays$RNA@counts)[1])
# pbmc.markers$qvalue <- qvalue
# top <- pbmc.markers[pbmc.markers$qvalue<0.05,]
# top10 <- top %>% top_n(n = 50, wt = abs(avg_logFC))
# top10_genes<- top10$gene
# temp <- temp_object[,temp_object$Cluster%in%c(ident.1,ident.2)]
# temp$Cluster <- as.factor(as.vector(temp$Cluster))

# # debugonce(annotated_heat)
# annotated_heat(object=temp,
#                   row_annotation=c(1),
#                   gene_list=top10_genes,
#                   Rowv=TRUE,
#                   gene_list_name="DF_genes",
#                   title=title,
#                   ordering="Cluster")

# DefaultAssay(temp_object) <- "integrated"
# write.table(pbmc.markers, file = paste0(Sys.Date(),"_TO_EXP_each_",target,"_",title,".tsv"),row.names=FALSE, na="", sep="\t")

# test <- split(top$gene,top$cluster)

# x=compareCluster(test, fun='enrichGO', OrgDb='org.Hs.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05,ont="BP")
# pdf(paste0(Sys.Date(),"_enrichGO_BP_",target,"_",title,".pdf"),width=12,height=10)
# dotplot(x, showCategory=5, includeAll=FALSE)
# dev.off()

# x=compareCluster(test, fun='enrichGO', OrgDb='org.Hs.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05,ont="MF")
# pdf(paste0(Sys.Date(),"_enrichGO_MF_",target,"_",title,".pdf"),width=12,height=10)
# dotplot(x, showCategory=5, includeAll=FALSE)
# dev.off()
# pbmc.markers[c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3"),]
# setwd("../")




































# dir.create("Reanalysis")
# setwd("Reanalysis")
# #Remove Oligos from Microglia
# Micro_Sub2 <- object_subset(object=Micro_Sub,Conditions=NULL,Clusters=c(1,2,3),Cell_Types=c("Microglia"),re_project=FALSE)
# Micro_Sub2$Cluster<-as.factor(as.vector(Micro_Sub2$Cluster))
# Micro_Sub2 <- reduce_dim(Micro_Sub2,project=project)$Combined

# p1 <- plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=TRUE,ncol=1)
# p2 <-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=TRUE,ncol=1)
# p3 <-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1 )
# pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# pdf(paste(Sys.Date(),project,"Cell_Type_tsne.pdf",sep="_"))
# print(p3)
# dev.off()
# pdf(paste(Sys.Date(),project,"Condition_tsne.pdf",sep="_"))
# print(p1)
# dev.off()
# pdf(paste(Sys.Date(),project,"Cluster_tsne.pdf",sep="_"))
# print(p2)
# dev.off()

# dir.create("DF_Clusters")
# setwd("DF_Clusters")
# return_fun <- df_genes(Micro_Sub2,"Cluster",top_n=50,logfc.threshold=0.8 ,n_cores=4,latent.vars=c("nCount_RNA"))
# setwd("../")









# gene_list <- c("JUN","JUND","FOSB","PLP1")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster")

# cells_list <- c("AAGAATAGGATC_1",
# "GCGCTTTGTGGA_1",
# "TAAATCCGAGTC_1",
# "CCACACTTTACG_1",
# "CAAGCGAAGTGA_1",
# "GTACTTAATTCC_1",
# "TAGCCCCAGAGT_1",
# "AGCGTGGGGTAT_1",
# "CCTCCGAACTGA_1",
# "TTTACACTTAGG_2",
# "TAAAGACTCGCG_2",
# "GCCATGCGAGGG_2",
# "GGTGTTTGAATC_2",
# "ACTCAGGCACTA_2",
# "GCCATGCGAGGA_2",
# "GCCATGCGAGGC_2",
# "GCCATGCGAGGT_2")

# gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster",title="Homeostatic")

# counts_mat <- as.matrix(Micro_Sub@assays$RNA@counts)
# intergraded_mat <- as.matrix(Micro_Sub@assays$integrated@data)
# scaled_mat <- Micro_Sub@assays$integrated@scale.data
# clusters <- Micro_Sub$Cluster
# mat <- counts_mat
# mat <- intergraded_mat

# mat2 <- as.data.frame(cbind(t(mat[rownames(mat)%in%gene_list,]),as.factor(clusters)))

# plot_list <- list()
# for (gene in gene_list){
# 	p <- ggplot(mat2, aes(x=clusters , y=gene )) +
#   		geom_violin()
#   	plot_list[[gene]] <- p
# }
# grid.arrange(grobs =plot_list,nrow=6)

# p <- ggplot(mat2, aes(x=clusters , y=CSF1R )) +
#   geom_violin()

# p <- ggplot(mat2, aes(x=clusters , y=CSF1R )) +
#   geom_violin()

# p <- ggplot(mat2, aes(x=clusters , y=CSF1R )) +
#   geom_violin()



# # =================================== DF Conditions =================================================
# dir.create("DF_Conditions2")
# setwd("DF_Conditions2")
# DefaultAssay(Combined) <- "RNA"
# return_fun <- df_genes(Micro_Sub,"condition",top_n=15,logfc.threshold=0.15 ,n_cores=4)
# DefaultAssay(Combined) <- "integrated"
# setwd("../")
# # ------------------------------------------------------------------------------------------------



# # =================================== DF Clusters =================================================
# dir.create("DF_Clusters")
# setwd("DF_Clusters")
# return_fun <- df_genes(Micro_Sub,"Cluster",top_n=15,logfc.threshold=0.15 ,n_cores=4)
# setwd("../")



# library(stackoverflow)
# cl_tops <- chunk2(return_fun$top_markers,length(unique(Micro_Sub$Cluster)))
# list_bar_cl <- c()
# for (i in cl_tops){
# 	list_bar_cl <- c(list_bar_cl,i[1:3])
# }
# plot_bar_cells(list_bar_cl[1:3],object=Micro_Sub,save=TRUE,title="CL1",target="Cluster")
# plot_bar_cells(list_bar_cl[4:6],object=Micro_Sub,save=TRUE,title="CL2",target="Cluster")
# plot_bar_cells(list_bar_cl[6:9],object=Micro_Sub,save=TRUE,title="CL3",target="Cluster")




# dir.create("GSEA")
# setwd("GSEA")

# sub_degs <- return_fun$degs[return_fun$degs$p_val_adj<0.05 & return_fun$degs$avg_logFC>1.5,]
# sub_degs <- return_fun$degs[return_fun$degs$p_val_adj<0.05,]
# universe <- rownames(Combined@assays$RNA@counts)

# symbols = unlist(lapply(tolower(as.vector(sub_degs$gene)),simpleCap))
# test <- split(symbols,sub_degs$cluster)

# x=compareCluster(test, fun='enrichGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05,ont="BP")
# pdf("enrichGO.pdf",width=12,height=10)
# dotplot(x, showCategory=5, includeAll=FALSE)
# dev.off()

# x2=compareCluster(test, fun='groupGO', OrgDb='org.Mm.eg.db',keyType="SYMBOL",ont="BP",level=2)

# pdf("groupGO.pdf",width=12)
# dotplot(x2, showCategory=10, includeAll=FALSE)
# dev.off()

# list_symbols <- list()
# iter=0
# for (i in names(test)){
# 	iter=iter+1
# 	eg=bitr(as.vector(test[[i]]),fromType = "SYMBOL",toType=c("SYMBOL","ENTREZID"),OrgDb='org.Mm.eg.db')$ENTREZID
# 	list_symbols[[i]] <-c(eg)
# }

# x3=compareCluster(list_symbols, fun='enrichKEGG', organism="mmu",pAdjustMethod = "BH",
#                         pvalueCutoff  = 0.05)

# pdf("enrichKEGG.pdf",width=12)
# dotplot(x3, showCategory=10, includeAll=FALSE)
# dev.off()



# # debugonce(gene_set_enrich)
# sub_degs <- return_fun$degs[return_fun$degs$p_val_adj<0.05 & return_fun$degs$avg_logFC>1.5,]
# universe <- rownames(Combined@assays$RNA@counts)
# for (cl in unique(sub_degs$cluster)){
# 	dir.create(paste0("Cluster_",cl))
# 	setwd(paste0("Cluster_",cl))
# 	gene_symbols <- sub_degs[sub_degs$cluster==cl,]$gene
# 	ego <-gene_set_enrich(gene_symbols,universe=universe,organism,ontology=c("BP","MF"),qval_thres=0.05,save=TRUE,title=cl)
# 	# pdf(paste(Sys.Date(),"GSEA_Dotplot","Cluster",cl,".pdf",sep="_"),height=14,width=18)
# 	# print(barplot(ego, showCategory=30))
# 	# dev.off()
# 	das <- david_enrich(gene_symbols=gene_symbols,
# 		idents=NULL,
# 		universe=universe,organism=organism,
# 		qval_thres=0.05,save=TRUE,title=paste0("CL",cl))
# 	setwd("../")
# }

# setwd("../")


# setwd("../")
# # ------------------------------------------------------------------------------------------------

# # =================================== Pseudotime =================================================
# dir.create("Pseudotime")
# setwd("Pseudotime")
# # debugonce(df_pseudotime)
# astro_pseudo <-df_pseudotime(Micro_Sub,reducedDim="DiffMap",reverse=FALSE,method="slingshot")
# setwd("../")
# # ------------------------------------------------------------------------------------------------
# setwd("../")
# # -----------------------------------------------------------------------------------------------


# find . -type f -name '*.pdf' -print0 |
#   while IFS= read -r -d '' file
#     do convert -verbose -density 300  "${file}" "TIFF/${file%.*}.tiff"
#   done

# convert ../*.pdf  -set filename:fname '%t_tn' +adjoin '%[filename:fname].tiff' -density=1200


# # ======================================== Hybrid Analysis  =================================
# dir.create("Hybrid_Analysis")
# setwd("Hybrid_Analysis")
# debugonce(object_subset)


# Hybrid_Sub <- object_subset(object=Combined,Conditions=NULL,Clusters=NULL,Cell_Types=c("Hybrid"))



# p1 <- plot_cells(Hybrid_Sub,target="condition",leg_pos="right",save=TRUE,ncol=1)
# p2 <-plot_cells(Hybrid_Sub,target="Cluster",leg_pos="right",save=TRUE,ncol=1)
# p3 <-plot_cells(Hybrid_Sub,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1 )
# pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# pdf(paste(Sys.Date(),project,"Cell_Type_tsne.pdf",sep="_"))
# print(p3)
# dev.off()
# pdf(paste(Sys.Date(),project,"Condition_tsne.pdf",sep="_"))
# print(p1)
# dev.off()
# pdf(paste(Sys.Date(),project,"Cluster_tsne.pdf",sep="_"))
# print(p2)
# dev.off()



# # =================================== DF Conditions =================================================
# dir.create("DF_Conditions")
# setwd("DF_Conditions")
# return_fun <- df_genes(Hybrid_Sub,"condition",top_n=15,logfc.threshold=0.25 ,n_cores=4)
# setwd("../")
# # ------------------------------------------------------------------------------------------------
# # =================================== DF Clusters =================================================
# dir.create("DF_Clusters")
# setwd("DF_Clusters")
# return_fun <- df_genes(Hybrid_Sub,"Cluster",top_n=15,logfc.threshold=0.25 ,n_cores=4)
# setwd("../")
# # ------------------------------------------------------------------------------------------------

# # =================================== Pseudotime =================================================
# dir.create("Pseudotime")
# setwd("Pseudotime")
# # debugonce(df_pseudotime)
# astro_pseudo <-df_pseudotime(Hybrid_Sub,reducedDim="DiffMap",reverse=FALSE,method="slingshot")
# setwd("../")
# # ------------------------------------------------------------------------------------------------
# setwd("../")
# # -----------------------------------------------------------------------------------------------















# # setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-11-20_seurat_elbow_FALSE")
# # Combined <- readRDS("Cell_Types_2019-11-20_seurat_elbow_FALSE.rds")
# # setwd("../")

# dir.create("Rescue")
# setwd("Rescue")
# Rescue_counts <- as.matrix(Combined@assays$RNA@counts)
# Corrected <- ComBat(Rescue_counts,Combined$condition)
# combat_edata2 = ComBat(dat=Rescue_counts, batch=Combined$condition, mod=NULL, par.prior=FALSE, mean.only=TRUE)
# f.out2 <- fastMNN(Rescue_counts, batch=Combined$condition)
# f.out2 <- logNormCounts(f.out2)
# manno <- runPCA(f.out2)
# manno.seurat <- as.Seurat(manno, counts = "counts", data = "logcounts")

# Rescue <- CreateSeuratObject(counts = combat_edata2, project = "all",
#                                  min.cells = 5,
#                                  min.features = 200)

# dim(Rescue)
# dim(Combined)
# Rescue$stim <- Combined$stim[colnames(Combined)%in%colnames(Rescue)]
# Rescue$condition <- Combined$condition[colnames(Combined)%in%colnames(Rescue)]

# # ====================== Normalization ====================
# Rescue <- NormalizeData(object = Rescue, normalization.method = "LogNormalize", scale.factor = 10000)
# # ----------------------------------------------


# # ====== Identification of highly variable features (feature selection)
# Rescue <- FindVariableFeatures(object = Rescue, selection.method = "vst", nfeatures = 2000)

# #Rescue <- ScaleData(object = Rescue, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)

# # Identify the 10 most highly variable genes
# top10 <- head(x = VariableFeatures(object = Rescue), 10)

# # plot variable features with and without labels
# plot1 <- VariableFeaturePlot(object = Rescue)
# plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
# # CombinePlots(plots = list(plot1, plot2))
# pdf(paste(Sys.Date(),"_top10_var_feat_",Rescue$stim,".pdf"))
# print(plot2)
# dev.off()


# dir.create("Clustering")
# setwd("Clustering")
# Rescue <- reduce_dim(Rescue,project=project)$Combined

# pdf(paste(Sys.Date(),project,"tsne","Conditions.pdf",sep="_"))
# plot_cells(Rescue,target="condition",leg_pos="right",save=TRUE,ncol=1)
# dev.off()
# pdf(paste(Sys.Date(),project,"tsne","Cluster.pdf",sep="_"))
# # debugonce(plot_cells)
# plot_cells(Rescue,target="Cluster",leg_pos="right",save=TRUE,ncol=1)
# dev.off()


# plot_nFeatures(Rescue,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)
# plot_tot_mRNA(Rescue,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE)






# Rescue$Cell_Type <- as.vector(Rescue$Cluster)
# Rescue$Cell_Type[Rescue$Cell_Type==4] <- "Microglia"
# # ======================================== Microglia Analysis  =================================
# dir.create("Microglia_Analysis")
# setwd("Microglia_Analysis")
# Micro_Sub <- object_subset(object=Rescue,Conditions=NULL,Clusters=c(4),Cell_Types="Microglia" )

# gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE)



# gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster",title="Homeostatic")


# # =================================== DF Clusters =================================================
# dir.create("DF_Clusters")
# setwd("DF_Clusters")
# # debugonce(df_genes)
# return_fun <- df_genes(Micro_Sub,"Cluster",top_n=80,logfc.threshold=0.15 ,n_cores=4)
# setwd("../")




# dir.create("DF_Conditions")
# setwd("DF_Conditions")
# # debugonce(df_genes)
# return_fun <- df_genes(Micro_Sub,"condition",top_n=50,logfc.threshold=0.15 ,n_cores=4)
# setwd("../")








# # ===================== STEFANO
# setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-11-25_seurat_elbow_FALSE_remove_ribsomal_TRUE")
# #saveRDS(Oihane,"Cell_Types_2019-11-25_seurat_elbow_FALSE_remove_ribsomal_TRUE.rds")
# Oihane <- readRDS("Cell_Types_2019-11-25_seurat_elbow_FALSE_remove_ribsomal_TRUE.rds")
# project ="MBSYN"
# dataset <- project

# plot_cells(Oihane,target="Cell_Type")


# dir.create("Microglia")
# setwd("Microglia")
# # ===== Subsetting =====
# Micro_Sub <- object_subset(object=Oihane,Conditions=NULL,Clusters=NULL,Cell_Types=c("Microglia") )
# # ===== DM Reduction ===
# Micro_Sub <- reduce_dim(Micro_Sub,project=project,resolution=c(1.2))$Combined

# gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
# plot_bar_cells(gene_list,object=Micro_Sub,save=TRUE,target="Cluster")


# p1<-plot_cells(Micro_Sub,target="condition",leg_pos="right",save=TRUE,ncol=1,reduction="umap")
# p2<-plot_cells(Micro_Sub,target="Cluster",leg_pos="right",save=TRUE,ncol=1,reduction="umap")
# p3<-plot_cells(Micro_Sub,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,reduction="umap")

# pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()

# p1<-plot_cells(Micro_Sub,target="condition",leg_pos="right",save=TRUE,ncol=1,reduction="tsne")
# p2<-plot_cells(Micro_Sub,target="Cluster",leg_pos="right",save=TRUE,ncol=1,reduction="tsne")
# p3<-plot_cells(Micro_Sub,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,reduction="tsne")

# pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
# print(plot_grid(p1, p2,p3,ncol = 3))
# dev.off()




# plot_cells(Micro_Sub,target="Cluster")
# dir.create("DF_Clusters")
# setwd("DF_Clusters")



# # =============================== PAIRWISE DF ===============================================
# dir.create("DF_Pairwise")
# setwd("DF_Pairwise")
# Idents(Micro_Sub) <- Micro_Sub$Cluster
# cl_combinations <- combn(levels(Micro_Sub$Cluster),2)

# for(comb in 1:dim(cl_combinations)[2]){
#     idents <- as.numeric(cl_combinations[,comb])
#     ident.1 <- idents[1]
#     ident.2 <- idents[2]
#     DefaultAssay(Micro_Sub) <- "RNA"
#     Micro_Sub <- ScaleData(Micro_Sub)
#     pbmc.markers <- FindMarkers(object = Micro_Sub,ident.1=ident.1,ident.2=ident.2,
#                                            assay ="RNA",
#                                            logfc.threshold=0,
#                                            only.pos = TRUE,
#                                            test.use = "MAST",latent.vars = c("nCount_RNA"))
#     pbmc.markers$gene <- rownames(pbmc.markers)
#     top10 <- pbmc.markers %>% top_n(n = 50, wt = abs(avg_logFC))
#     top10_genes<- top10$gene
#     temp <- Micro_Sub[,Micro_Sub$Cluster%in%c(ident.1,ident.2)]
#     temp$Cluster <- as.factor(as.vector(temp$Cluster))
#     # debugonce(annotated_heat)
#     annotated_heat(object=temp,
#                       row_annotation=c(1),
#                       gene_list=top10_genes,
#                       gene_list_name="DF_genes",
#                       title=paste(cl_combinations[,comb],collapse = "_"),
#                       ordering="Cluster")

#     DefaultAssay(Micro_Sub) <- "integrated"

#     pbmc.markers[c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3"),]
# }




# Homeostatic <- c(grep("FCRLS",return_fun$degs$gene),
# grep("LGMN",return_fun$degs$gene),
# grep("MARCKS",return_fun$degs$gene),
# grep("CSF1R",return_fun$degs$gene),
# grep("P2RY12",return_fun$degs$gene),
# grep("SERINC3",return_fun$degs$gene))

# return_fun$degs[Homeostatic,]
# # ----------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------
pdf('DF_VLN.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = return_fun_cl_micro2$top_markers,group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()



setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-12-11_seurat_elbow_FALSE_remove_ribsomal_TRUE/Microglia_Analysis/Reanalysis_Counts/DF_Clusters")
pdf('DF_VLN_1.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = top$gene_short_name [1:9],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()

pdf('DF_VLN_2.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = top$gene_short_name[10:19],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()

pdf('DF_VLN_3.pdf',height=24)
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = top$gene_short_name[20:27],group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()

setwd("../")

pdf('VLN_APOE.pdf')
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = c("APOE"),group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
dev.off()


Homeostatic <- c(grep("FCRLS",pbmc.markers$gene),
grep("LGMN",pbmc.markers$gene),
grep("MARCKS",pbmc.markers$gene),
grep("CSF1R",pbmc.markers$gene),
grep("P2RY12",pbmc.markers$gene),
grep("SERINC3",pbmc.markers$gene))



 pbmc.markers[Homeostatic,]


Homeostatic <- c(grep("FCRLS",pbmc.markers2$gene),
grep("LGMN",pbmc.markers2$gene),
grep("MARCKS",pbmc.markers2$gene),
grep("CSF1R",pbmc.markers2$gene),
grep("P2RY12",pbmc.markers2$gene),
grep("SERINC3",pbmc.markers2$gene))

pbmc.markers2[Homeostatic,]




gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","FOSB","SERINC3")
plot_bar_cells(gene_list,object=Micro_Sub2,save=FALSE,target="Cluster",title="else")



dir.create("Paper_Figures")
setwd("Paper_Figures")
color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_cells <- c(rgb(255,0,0, maxColorValue = 255),rgb(0,0,255, maxColorValue = 255),rgb(0,0,160, maxColorValue = 255),rgb(255,100,177, maxColorValue = 255),rgb(128,0,128, maxColorValue = 255),rgb(128,64,0, maxColorValue = 255),rgb(0,255,0, maxColorValue = 255))
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)




p1 <- plot_cells(Combined,target="condition",leg_pos="right",save=TRUE,ncol=1,color_list=color_list)
p2 <-plot_cells(Combined,target="Cluster",leg_pos="right",save=TRUE,ncol=1,color_list=color_list)
p3 <-plot_cells(Combined,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,color_list=color_list)

pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

pdf(paste(Sys.Date(),project,"Cell_Type_tsne.pdf",sep="_"))
print(p3)
dev.off()
pdf(paste(Sys.Date(),project,"Condition_tsne.pdf",sep="_"))
print(p1)
dev.off()
pdf(paste(Sys.Date(),project,"Cluster_tsne.pdf",sep="_"))
print(p2)
dev.off()


p1 <- plot_cells(Combined,target="condition",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list)
p2 <-plot_cells(Combined,target="Cluster",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list)
p3 <-plot_cells(Combined,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list)

pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

pdf(paste(Sys.Date(),project,"Cell_Type_umap.pdf",sep="_"))
print(p3)
dev.off()
pdf(paste(Sys.Date(),project,"Condition_umap.pdf",sep="_"))
print(p1)
dev.off()
pdf(paste(Sys.Date(),project,"Cluster_umap.pdf",sep="_"))
print(p2)
dev.off()

color_clust_microglia <- c(rgb(128, 128, 192, maxColorValue = 255),
    rgb(255,128,192, maxColorValue = 255),
    rgb(128,0,64, maxColorValue = 255),
    rgb(128,0,255, maxColorValue = 255))
color_cells_microglia <- c(rgb(128,0,128, maxColorValue = 255))
color_list_microglia <- list(condition=color_cond,Cluster=color_clust_microglia,Cell_Type=color_cells_microglia,State=color_clust_microglia)




dir.create("Micro_Sub2")
setwd("Micro_Sub2")

p1 <- plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=TRUE,ncol=1,color_list=color_list_microglia)
p2 <-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=TRUE,ncol=1,color_list=color_list_microglia)
p3 <-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,color_list=color_list_microglia)

pdf(paste(Sys.Date(),project,"tsne.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

pdf(paste(Sys.Date(),project,"Cell_Type_tsne.pdf",sep="_"))
print(p3)
dev.off()
pdf(paste(Sys.Date(),project,"Condition_tsne.pdf",sep="_"))
print(p1)
dev.off()
pdf(paste(Sys.Date(),project,"Cluster_tsne.pdf",sep="_"))
print(p2)
dev.off()

p1 <- plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list_microglia)
p2 <-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list_microglia)
p3 <-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=TRUE,ncol=1,reduction="umap",color_list=color_list_microglia)

pdf(paste(Sys.Date(),project,"umap.pdf",sep="_"),width=16,height=8)
print(plot_grid(p1, p2,p3,ncol = 3))
dev.off()

pdf(paste(Sys.Date(),project,"Cell_Type_umap.pdf",sep="_"))
print(p3)
dev.off()
pdf(paste(Sys.Date(),project,"Condition_umap.pdf",sep="_"))
print(p1)
dev.off()
pdf(paste(Sys.Date(),project,"Cluster_umap.pdf",sep="_"))
print(p2)
dev.off()


pdf('DF_Microsub2_Violin.test_RNA.pdf')
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = c("FOSB","JUNB","LGMN"),group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()

pdf('APOE_Microsub2_Violin.test_RNA.pdf')
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = c("APOE"),group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
CombinePlots(plots = plots, legend = 'right',ncol =1 )
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()

setwd("../")

setwd("../")




df_qc <- as.data.frame(Micro_Sub2@assays$RNA@counts[c("P2RY12","LGMN"),],Cluster=Micro_Sub2$Cluster)
p1 <- ggplot(,aes(factor(Micro_Sub2$Cluster),df_qc[1,],fill=Micro_Sub2$Cluster))+
    geom_boxplot(width=0.8,show.legend=FALSE)+
    xlab("")+
    ylab("Expression Level")





# =============== 2020 ==========================================
setwd("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-12-11_seurat_elbow_FALSE_remove_ribsomal_TRUE")
load("Paper_WORKSPACE.RData")

dir.create("2020_oihane_extra_plots")
setwd("2020_oihane_extra_plots")

Micro_Sub2 <- Micro_Sub
color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
color_clust_microglia <- c("#42858C",rgb(128, 128, 192, maxColorValue = 255),
                           rgb(128,0,255, maxColorValue = 255),
                           rgb(255,128,192, maxColorValue = 255),
                           rgb(128,0,64, maxColorValue = 255))

p1 <- DimPlot(Micro_Sub2,group.by = "Cluster",cols= color_clust_microglia )

Micro_Sub2 <- object_subset(object=Micro_Sub,Conditions=NULL,Clusters=c(1,2,3,4,6,7),Cell_Types=NULL)
DefaultAssay(Micro_Sub2) <- "RNA"
Micro_Sub2 <- reduce_dim(Micro_Sub2,project=project,resolution=c(1.1))$Combined
DimPlot(Micro_Sub2,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)


temp <- as.numeric(Micro_Sub2$Cluster)
temp[temp==4] <- 5
temp[temp==3] <- 4
temp[temp==5] <- 3
Micro_Sub2$Cluster <- as.factor(temp)
DimPlot(Micro_Sub2,group.by = "Cluster")
Micro_Sub2$micsub <- as.vector(Micro_Sub2$Cluster)

Micro_Sub2$micsub[Micro_Sub2$micsub==1]<- "1.Homeostatic subset"
Micro_Sub2$micsub[Micro_Sub2$micsub==2]<- "2.Intermediate subset 1"
Micro_Sub2$micsub[Micro_Sub2$micsub==3]<- "3.Intermediate subset 2"
Micro_Sub2$micsub[Micro_Sub2$micsub==4]<- "4.Immune alerted subset"
levels(Micro_Sub2$micsub)<- c("Homeostatic subset", "Intermediate subset 1", "Intermediate subset 2", "Immune alerted subset")


DefaultAssay(Micro_Sub2) <- "RNA"
new_names <- unlist(lapply(tolower(rownames(Micro_Sub2@assays$RNA@data)),simpleCap))
rownames(Micro_Sub2@assays$RNA@data) <- new_names
features <- c("SOCS3","CD14","GPR84","LYZ1","CASP4","FTH1",
              "ICAM1","IL1B","ADAMTS1","CD83","CCL4","NFKBIZ")
s_features <- factor(unlist(lapply(tolower(features),simpleCap)))
p3<-DotPlot(object = Micro_Sub2, features = s_features,group.by = "micsub") +
    RotatedAxis()+
    scale_x_discrete(limits=s_features)

features <- c("H2.AA","H2.AB1","CD74")
s_features <- unlist(lapply(tolower(features),simpleCap))
p4<-DotPlot(object = Micro_Sub2, features = s_features,group.by = "micsub") +
    RotatedAxis()+
    scale_x_discrete(limits=s_features)

features <- c("HEXB","CX3CR1","P2RY12","C1QA","FCRLS")
s_features <- unlist(lapply(tolower(features),simpleCap))
p5<-DotPlot(object = Micro_Sub2, features = s_features,group.by = "micsub") +
    RotatedAxis()+
    scale_x_discrete(limits=s_features)
#DotPlot(object = Micro_Sub2, features = s_features,group.by = "micsub")+ theme(axis.text.x = element_text(angle = 45, hjust=1),axis.text.y = element_text(angle =30, hjust=1))

p3

############

library(viridis)

pdf("Scatter.pdf")
p1
dev.off()

pdf("Heatmap.pdf",width=12)
DoHeatmap(Micro_Sub2,features = degs_monocle_plot,group.by = "Cluster", label=FALSE,
          group.colors = color_clust_microglia,size = 0.3)+
    scale_fill_gradientn(colors = rev(brewer.pal(name = "RdYlBu",n =110)))+
    theme(text = element_text(size = 10))
dev.off()


pdf("DotPlots1.pdf",width=12)
p3
dev.off()

pdf("DotPlots2.pdf",width=8)
p4
dev.off()

pdf("DotPlots3.pdf",width=10)
p5
dev.off()



#' Annotated Heatmap
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: object class from Monocle/Seurat.
#' @param row_annotation: a vector with gene marker category
#' @param gene_list: A vector of gene names that want to score.
#' @param ordering:  Preferable order of cells (condition/Cluster)
#' @param Colv : Hierarchical clustering of columns (Default = NA)
#' @param Rowv : Hierarchical clustering of rows (Default = NA)
#'
#' @return An annotated heatmap
annotated_heat<- function(object,row_annotation,gene_list,gene_list_name,title,ordering="Cluster",Colv=NA,Rowv=NA,One_annot=FALSE){
    tool <- object_identifier(object)
    if (toupper(ordering)=="CLUSTER"){
        temp_ordering  <- order(object$Cluster,object$condition)
    } else{
        temp_ordering  <- order(object$condition,object$Cluster)
    }
    Condi <- object$condition[temp_ordering]
    clusts    <- object$Cluster[temp_ordering]
    r_annotation = row_annotation
    if (ordering=="Cluster"){
        c_annotation = data.frame(Clusters = clusts,Conditions=Condi)
    } else {
        c_annotation = data.frame(Conditions=Condi,Clusters = clusts)
    }

    if (One_annot){
        c_annotation2 <- data.frame(c_annotation[,1])
        colnames(c_annotation2) <- colnames(c_annotation)[1]
        c_annotation <- c_annotation2
    }



    if(tolower(tool=="monocle")){
        object_rownames <- fData(object)$gene_short_name
        GE_matrix<- data.matrix(object, rownames.force = NA)
        Normilize_Vector <- pData(object)[, 'Size_Factor']
        Norm1 <- as.matrix(sweep(GE_matrix,2,Normilize_Vector,"/"))
        GE_matrix <- as.matrix(Norm1)

    }else{
        GE_matrix <-as.matrix(object@assays[[object@active.assay]]@data)
        GE_matrix<- GetAssayData(object = object, slot = "scale.data")
        GE_matrix <- as.matrix(GE_matrix)
        object_rownames <-rownames(GE_matrix)

    }
    disp.min <- -2.5
    disp.max <- 2.5
    test0 <- MinMax(data = GE_matrix, min = disp.min, max = disp.max)
    #test0 <- GE_matrix

    list_ind1 <- match(gene_list,object_rownames)

    na_index <- is.na(list_ind1)
    r_annotation<-data.frame(gene_list_name=r_annotation[!na_index])
    gene_list <- gene_list[!na_index]
    list_ind1<- list_ind1[!na_index]

    # Subset Top genes
    Top_DF_Cond <- test0[list_ind1,]
    rownames(Top_DF_Cond)<-gene_list



    color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
    color_clust <- c("#42858C",rgb(128, 128, 192, maxColorValue = 255),
                               rgb(128,0,255, maxColorValue = 255),
                               rgb(255,128,192, maxColorValue = 255),
                               rgb(128,0,64, maxColorValue = 255))
    color_cells <- c(rgb(255,0,0, maxColorValue = 255),rgb(0,0,255, maxColorValue = 255),rgb(0,0,160, maxColorValue = 255),rgb(255,100,177, maxColorValue = 255),rgb(128,0,128, maxColorValue = 255),rgb(128,64,0, maxColorValue = 255),rgb(0,255,0, maxColorValue = 255))
    color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

    colr_annot <- list(gene_list_name=color_cells[1:length(unique(row_annotation))],
                       Conditions=color_cond[1:length(levels(Condi))],
                       Clusters=color_clust[1:length(levels(clusts))])

    if(length(unique(row_annotation))==1){
        aheatmap(as.matrix(Top_DF_Cond)[,temp_ordering],
                 annCol = c_annotation,
                 Colv = Colv,
                 Rowv = Rowv,
                 annColors=colr_annot,
                 labRow = NULL,
                 #scale="r1",
                 # color=CustomPalette(low = "#000000", high = "#ff0000", k = 100),
                 cexRow=3.5,
                 color=rev(brewer.pal(name = "RdYlBu",n =11)),#CustomPalette(low = "#ff0905", mid = "#000000" , high = "#ffff05", k = 100),
                 #color=CustomPalette(low = "#000000" , high = "#db1200", mid = "#400459", k = 50),
                 #PurpleAndYellow(),#rev(brewer.pal(name = "RdYlBu",n =11)),
                 #color=jcolors(palette = "pal11"),
                 main=paste("Gene marker for ",title),
                 filename=paste(Sys.Date(),"_Heat_",title,"_target-",ordering,".pdf",sep=""),width = 13, height = 10)
    }else{
        aheatmap(as.matrix(Top_DF_Cond)[,temp_ordering],
                 annCol = c_annotation,
                 annRow= r_annotation,
                 Colv = Colv,
                 Rowv = Rowv,
                 annColors=colr_annot,
                 cexRow=3.5,
                 color=rev(brewer.pal(name = "RdYlBu",n =11)),#rev(brewer.pal(name = "RdYlBu",n =11)),#,jcolors(palette = "pal12"),
                 labRow = NULL,
                 main=paste("Gene marker for ",title),
                 filename=paste(Sys.Date(),"_Heat_",title,"_target-",ordering,".pdf",sep=""),width = 13, height = 10)
    }

}



annotated_heat(Micro_Sub2,c(1),degs_monocle_plot,"Cluster","Cluster",
               ordering="Cluster",Colv=NA,Rowv=NA,One_annot=FALSE)



library(ggpubr)
ggarrange(ggarrange(p1,p2),ggarrange(p3,p4,p5,nrow=1),ncol=1)



gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","MEG3","CALD1")
plot_bar_cells(gene_list,object=Combined,save=TRUE,title="all",normalized=FALSE,color_list=color_list)
plot_bar_cells(gene_list,object=Combined,save=TRUE,title="all",normalized=TRUE,color_list=color_list)






gene_list <- c("APOE","TREM2","TLR2","REEP4","SPP1","CTSD","CTSB","TLR4")
plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster",normalized=TRUE)
plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster",normalized=FALSE)
boxplot_ics(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list,normalized=TRUE)
boxplot_ics(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list,normalized=FALSE)
jitter_ics(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list,normalized=TRUE)
jitter_ics(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list,normalized=FALSE)
violin_ics(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list,normalized=TRUE)
violin_ics(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list,normalized=FALSE)
VlnPlot(object=Micro_Sub2,target="Cluster",features=gene_list,color_list=color_list)



gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","FOSB","SERINC3")
plot_bar_cells(gene_list,object=Micro_Sub2,save=FALSE,target="Cluster",title="microsub2",normalized=FALSE)
plot_bar_cells(gene_list,object=Micro_Sub2,save=FALSE,target="Cluster",title="microsub2",normalized=TRUE)
RidgePlot(object = Micro_Sub2, features = gene_list)
DotPlot(object = Micro_Sub2, features = gene_list,cols=brewer.pal(9, "OrRd"))


pdf("Gene_exp.pdf")
for (gene in gene_list){
    print(FeaturePlot(Micro_Sub2, features = gene,
            reduction = "tsne",
            pt.size=1.2,,cols=brewer.pal(9, "OrRd")))
}
dev.off()



data <- data.frame(cbind(factor(Micro_Sub2$Cluster),paste(" ",Micro_Sub2$condition)))
colnames(data)<- c("Cluster","Condition")
pdf(paste(Sys.Date(),'Microglia2_Barplot_num-Cond_per_Cluster.pdf',sep="_"))
print(ggplot(data, aes(Cluster)) +  geom_bar(aes(fill = Condition))+
          guides(col=guide_legend(ncol=2,))+
          scale_fill_manual(values =color_list_microglia$condition)+
          theme(legend.position="bottom") )
dev.off()



data <- data.frame(cbind(Combined$Cell_Type,paste(" ",Combined$condition)))

colnames(data)<- c("Cell_Type","Condition")



pdf(paste(Sys.Date(),'Cell_Type_Barplot_num-Cond_per_Cluster.pdf',sep="_"))
print(ggplot(data, aes(Cell_Type)) +  geom_bar(aes(fill = Condition))+
          guides(col=guide_legend(ncol=2,))+
          scale_fill_manual(values =color_list$condition)+
          theme(legend.position="bottom") )

dev.off()


dir.create("DF_Conditions")
setwd("DF_Conditions")
return_fun_cond_micro2 <- df_genes(Micro_Sub2,"condition",top_n=15,logfc.threshold=0,
    n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE,color_list=color_list_microglia)

micro2_all_cond_df <- return_fun_cond_micro2$degs$gene[return_fun_cond_micro2$degs$p_val_adj<0.05]
annotated_heat(Micro_Sub2,c(1),micro2_all_cond_df,"condition","condition",
    ordering="condition",Colv=NA,Rowv=NA,One_annot=FALSE,color_list=color_list_microglia)

setwd("../")
# -------------

# ======= Could you generate t-SNE graphs with the expression of the 78 DEGs between microglia subsets?
# Using Monocle DFs_Count Clustering

file <- "C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-12-11_seurat_elbow_FALSE_remove_ribsomal_TRUE/microglia-heterogeneity/Microglia_Analysis/Reanalysis_Counts/DF_Clusters/monocle"

degs_monocle_plot <- read.table(paste0(file,"/Monocle_Top_qvalue.tsv"))
degs_monocle_plot<- as.vector(degs_monocle_plot[degs_monocle_plot$qval<0.05,"gene_short_name"])



pdf('Scatter_DF_Micro2_Mon_data_umap.pdf')
plots <-FeaturePlot(object = Micro_Sub2,features = degs_monocle_plot,ncol = 1,combine=FALSE,reduction="umap",
            pt.size=1.2,cols=brewer.pal(9, "OrRd")[-1])
plots <- lapply(
  X = plots,
  FUN = function(p) print(p)
)
dev.off()

pdf('Scatter_DF_Micro2_Mon_data_tsne.pdf')
plots <-FeaturePlot(object = Micro_Sub2,features = degs_monocle_plot,ncol = 1,combine=FALSE,reduction="tsne",
            pt.size=1.2,cols=brewer.pal(9, "OrRd")[-1])
plots <- lapply(
  X = plots,
  FUN = function(p) print(p)
)
dev.off()


pdf('Vln_DF_Micro2_Mon_data_tsne.pdf')
plots <-VlnPlot(object = Micro_Sub2,features = degs_monocle_plot,ncol = 1,combine=FALSE,group.by="Cluster",
            cols=color_list_microglia$Cluster)
plots <- lapply(
  X = plots,
  FUN = function(p) print(p)
)
dev.off()

pdf('dOT_DF_Micro2_Mon.pdf')
# DotPlot(object = Micro_Sub2,features = degs_monocle_plot[1:35],group.by="Cluster",cols=brewer.pal(9, "OrRd"))
DotPlot(object = Micro_Sub2,features = degs_monocle_plot[1:20],group.by="Cluster")+ RotatedAxis()
DotPlot(object = Micro_Sub2,features = degs_monocle_plot[21:40],group.by="Cluster")+ RotatedAxis()
DotPlot(object = Micro_Sub2,features = degs_monocle_plot[41:60],group.by="Cluster")+ RotatedAxis()
DotPlot(object = Micro_Sub2,features = degs_monocle_plot[61:78],group.by="Cluster")+ RotatedAxis()
dev.off()

RidgePlot(pbmc, features = features, ncol = 2)

pdf('Ridge_DF_Micro2_Mon.pdf')
# DotPlot(object = Micro_Sub2,features = degs_monocle_plot[1:35],group.by="Cluster",cols=brewer.pal(9, "OrRd"))
RidgePlot(object = Micro_Sub2,features = degs_monocle_plot[1:20],group.by="Cluster")+ RotatedAxis()
RidgePlot(object = Micro_Sub2,features = degs_monocle_plot[21:40],group.by="Cluster")+ RotatedAxis()
RidgePlot(object = Micro_Sub2,features = degs_monocle_plot[41:60],group.by="Cluster")+ RotatedAxis()
RidgePlot(object = Micro_Sub2,features = degs_monocle_plot[61:78],group.by="Cluster")+ RotatedAxis()
dev.off()

annotated_heat(Micro_Sub2,c(1),degs_monocle_plot,"Cluster","Cluster",
    ordering="Cluster",Colv=NA,Rowv=NA,One_annot=FALSE,color_list=color_list_microglia)





# pdf('Scatter_DF_Micro2_Mon_counts.pdf')
# for(gene_name in degs_monocle_plot){
#     gene <- FetchData(Micro_Sub2,vars=gene_name,slot = "counts")[,1]
#     df <- as.data.frame(cbind(Micro_Sub2@reductions$umap@cell.embeddings,Counts=gene))
#     b <- ggplot(df, aes(x = UMAP_1, y = UMAP_2))
#     b <- b + geom_point(aes(color = Counts), size = 3) +
#       ggtitle(simpleCap(tolower(gene_name)))+
#       scale_color_gradientn(colors =  brewer.pal(9, "OrRd")) +
#       theme_cowplot()+theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
#     print(b)
# }
# dev.off()


# pdf('Scatter_DF_Micro2_Mon_data.pdf')
# for(gene_name in degs_monocle_plot){
#     gene <- FetchData(Micro_Sub2,vars=gene_name,slot = "data")[,1]
#     df <- as.data.frame(cbind(Micro_Sub2@reductions$umap@cell.embeddings,Counts=gene))
#     b <- ggplot(df, aes(x = UMAP_1, y = UMAP_2))
#     b <- b + geom_point(aes(color = Counts), size = 3) +
#       ggtitle(simpleCap(tolower(gene_name)))+
#       scale_color_gradientn(colors =  brewer.pal(9, "OrRd")) +
#       theme_cowplot()+theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
#     print(b)
# }
# dev.off()

# FeaturePlot(object = Micro_Sub2,features = gene_name,reduction="umap",slot="counts",
#             pt.size=1.2,cols=brewer.pal(9, "OrRd")[-1])


dir.create("Papaer_Figures_Numbered")
setwd("Papaer_Figures_Numbered")

cell_types <- Combined$Cell_Type
cell_types <- as.vector(Combined$Cluster)
save_cell_types<-Combined$Cell_Type
cell_types[ which(Combined$Cluster %in% c(1))] <- "Astrocytes"
cell_types[ which(Combined$Cluster %in% c(2))] <- "Microglia"
cell_types[ which(Combined$Cluster %in% c(3))] <- "Oligodendrocytes"
cell_types[ which(Combined$Cluster %in% c(4))] <- "Endothelial"
cell_types[ which(Combined$Cluster %in% c(5))] <- "Hybrid 1"
cell_types[ which(Combined$Cluster %in% c(6))] <- "Ependymal"
cell_types[ which(Combined$Cluster %in% c(7))] <- "Epithelial"
cell_types[ which(Combined$Cluster %in% c(8))] <- "Neurons"
cell_types[ which(Combined$Cluster %in% c(9))] <- "Hybrid 2"
Combined$Cell_Type <-cell_types

color_clust_ord_cl <- color_clust[c(1,4,6,7,5,9,2,8,3)]
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_clust_ord_cl,State=color_clust)


# ======================== PLOT 1 ============================================
pdf(paste("1",Sys.Date(),project,"Data_projections_BCE.pdf",sep="_"))
p<-plot_cells(Combined,target="condition",leg_pos="right",save=FALSE,ncol=1,color_list=color_list)
p<-plot_cells(Combined,target="Cluster",leg_pos="right",save=FALSE,ncol=1,color_list=color_list)
p<-plot_cells(Combined,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,color_list=color_list)
p<-plot_cells(Combined,target="condition",reduction="umap",leg_pos="right",save=FALSE,ncol=1,color_list=color_list)
p<-plot_cells(Combined,target="Cluster",reduction="umap",leg_pos="right",save=FALSE,ncol=1,color_list=color_list)
p<-plot_cells(Combined,target="Cell_Type",reduction="umap",leg_pos="right",save=FALSE,ncol=1,color_list=color_list)
dev.off()



pdf(paste("2",Sys.Date(),project,"Scatter_Markers_FGHI.pdf",sep="_"))
DefaultAssay(Combined) <- "RNA"
Combined <- ScaleData(Combined,, features = rownames(Combined))
FeaturePlot(Combined,features=c("PLPP3", "SLC1A2", "GJA1", "AQP4"),reduction="tsne",pt.size=2,cols=brewer.pal(9, "OrRd"))
FeaturePlot(Combined,features=c("HEXB", "P2RY12", "CX3CR1", "SIGLECH"),reduction="tsne",pt.size=2,cols=brewer.pal(9, "OrRd"))
FeaturePlot(Combined,features=c( "PLP1", "MBP", "MOBP", "TRF"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
FeaturePlot(Combined,features=c("LY6C1", "CLDN5", "PLTP", "PECAM1"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
FeaturePlot(Combined,features=c("SCN7A", "MAP2", "C1QL1", "PDGFRA"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
FeaturePlot(Combined,features=c("CCDC153", "TMEM212", "DYNLRB2", "RSPH4A"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
FeaturePlot(Combined,features=c("KCNJ13", "FOLR1", "CLIC6", "KL"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
FeaturePlot(Combined,features=c( "MEG3", "SNHG11", "NDRG4", "ENO2"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
FeaturePlot(Combined,features=c( "CALD1", "VTN", "NOTCH3", "SNAP25"),reduction="tsne",cols=brewer.pal(9, "OrRd"),pt.size=2)
dev.off()

# ================ 3
pie_per_CellType(Combined)
plot_bar_cells(Combined)
# ============================= MICROGLIA SUBPOPULATION COUNTS REPROJECTION ==================================================
pdf(paste("4",Sys.Date(),project,"Microglia2_projections.pdf",sep="_"))
p<-plot_cells(Micro_Sub2,target="condition",leg_pos="right",save=FALSE,ncol=1,color_list=color_list_microglia,pt.size=2)
p<-plot_cells(Micro_Sub2,target="Cluster",leg_pos="right",save=FALSE,ncol=1,color_list=color_list_microglia,pt.size=2)
p<-plot_cells(Micro_Sub2,target="Cell_Type",leg_pos="right",save=FALSE,ncol=1,color_list=color_list_microglia,pt.size=2)
p<-plot_cells(Micro_Sub2,target="condition",reduction="umap",leg_pos="right",save=FALSE,ncol=1,color_list=color_list_microglia,pt.size=2)
p<-plot_cells(Micro_Sub2,target="Cluster",reduction="umap",leg_pos="right",save=FALSE,ncol=1,color_list=color_list_microglia,pt.size=2)
p<-plot_cells(Micro_Sub2,target="Cell_Type",reduction="umap",leg_pos="right",save=FALSE,ncol=1,color_list=color_list_microglia,pt.size=2)
dev.off()
# -------------------------------------------------------------------------------------------------

pdf(paste("5",Sys.Date(),project,"Microglia2_Inflammation-activation_markers.pdf",sep="_"))
gene_list <- c("NFKBIZ", "CCL4" ,"CD83", "ADAMTS1", "IL1B", "NFKBIA",  "IL1A", "ICAM1", "FTH1",  "CASP4", "LYZ1", "GPR84", "CD14","SOCS3") #"ITM2B",
DotPlot(object = Micro_Sub2,group.by="Cluster", features = gene_list)+ RotatedAxis()
gene_list <- c("FCRLS", "C1QA", "P2RY12", "CX3CR1", "HEXB")
DotPlot(object = Micro_Sub2,group.by="Cluster", features = gene_list)+ RotatedAxis()
gene_list <- c("CD74","H2.AB1", "H2.AA")
DotPlot(object = Micro_Sub2,group.by="Cluster", features = gene_list)+ RotatedAxis()
gene_list <- c("ATF3", "EGR1", "H3F3B", "BTG2", "GPX1",
 "CHD4", "ATP5A1", "SDC4", "VPS29", "KLF2")
DotPlot(object = Micro_Sub2,group.by="Cluster", features = gene_list)+ RotatedAxis()

gene_list <- c("JUNB", "FOSB", "FOS", "JUN", "JUND")
DotPlot(object = Micro_Sub2,group.by="Cluster", features = gene_list)+ RotatedAxis()
dev.off()

#. Scatter plots UMAP (bigger size of the dots):
pdf(paste("6",Sys.Date(),project,"Microglia2_scatter_markers.pdf",sep="_"))
gene_list <- c("CD74","H2.AB1", "H2.AA","H2.K1")
FeaturePlot(Micro_Sub2,features=gene_list,reduction="umap",cols=brewer.pal(9, "OrRd"),pt.size=2,order=T)
dev.off()

# ======================== PLOT 2 ============================================

# ====================== BARPLOT CELL BELONGS TO WHICH CLUSTERS ====================================
data <- data.frame(cbind(factor(Micro_Sub2$Cluster),paste(" ",Micro_Sub2$condition)))
colnames(data)<- c("Cluster","Condition")
pdf(paste("2",Sys.Date(),'Barplot_num-Cond_per_Cluster_Microglia2.pdf',sep="_"))
print(ggplot(data, aes(Cluster)) +  geom_bar(aes(fill = Condition))+
          guides(col=guide_legend(ncol=2,))+
          scale_fill_manual(values =color_list_microglia$condition)+
          theme(legend.position="bottom") )
dev.off()
setwd("../")
# -------------------------------------------------------------------------------------------------




# ======================== PLOT 3 ============================================

file <- "C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/MBSYN/2019-12-11_seurat_elbow_FALSE_remove_ribsomal_TRUE/microglia-heterogeneity/Microglia_Analysis/Reanalysis_Counts/DF_Clusters/monocle"

degs_monocle_plot <- read.table(paste0(file,"/Monocle_Top_qvalue.tsv"))
degs_monocle_plot<- as.vector(degs_monocle_plot[degs_monocle_plot$qval<0.05,"gene_short_name"])

annotated_heat(object=Micro_Sub2,row_annotation=c(1),gene_list=degs_monocle_plot,
    gene_list_name="Cluster",
    title="3_Cluster",ordering="Cluster",
    Colv=NA,Rowv=NA,One_annot=FALSE,
    color_list=color_list_microglia)

pdf(paste0("4.1_",Sys.Date(),'_Scatter_Micro2_Mon_DF_data_tsne.pdf'))
plots <-FeaturePlot(object = Micro_Sub2,features = degs_monocle_plot,ncol = 1,combine=FALSE,reduction="tsne",
            pt.size=1.2,cols=brewer.pal(9, "OrRd")[-1])
plots <- lapply(
  X = plots,
  FUN = function(p) print(p)
)
dev.off()
# -------------------------------------------------------------------------------------------------


# ======================== PLOT 4 ============================================

pdf(paste0("4.2_",Sys.Date(),'_Scatter_Micro2_Mon_DF_data_umap.pdf'))
plots <-FeaturePlot(object = Micro_Sub2,features = degs_monocle_plot,ncol = 1,combine=FALSE,reduction="umap",
            pt.size=1.2,cols=brewer.pal(9, "OrRd")[-1])
plots <- lapply(
  X = plots,
  FUN = function(p) print(p)
)
dev.off()
# -------------------------------------------------------------------------------------------------



# ======================== PLOT 5 ============================================

pdf('5.1_Vln_Micro2_Mon_DF_data.pdf')
for(gene_iter in degs_monocle_plot){
    print(jitter_ics(object = Micro_Sub2,target = "Cluster",features = gene_iter,type = "violin",color_list = color_list_microglia,organism="mouse"))
}
dev.off()

pdf('5.2_Box_Micro2_Mon_DF_data.pdf')
for(gene_iter in degs_monocle_plot){
    print(jitter_ics(object = Micro_Sub2,target = "Cluster",features = gene_iter,type = "box",color_list = color_list_microglia,organism="mouse"))
}
dev.off()
# -------------------------------------------------------------------------------------------------


# ======================== PLOT 6 ============================================
gene_list_apoe <- c("APOE","TREM2","TLR2","REEP4","SPP1","CTSD","CTSB","TLR4")

pdf('6.1_Vln_Micro2_Gene_Lists_data.pdf')
for(gene_iter in gene_list_apoe){
    print(jitter_ics(object = Micro_Sub2,target = "Cluster",features = gene_iter,type = "violin",color_list = color_list_microglia,organism="mouse"))
}
dev.off()

pdf('6.2_Box_Micro2_Gene_Lists_data.pdf')
for(gene_iter in gene_list_apoe){
    print(jitter_ics(object = Micro_Sub2,target = "Cluster",features = gene_iter,type = "box",color_list = color_list_microglia,organism="mouse"))
}
dev.off()

# -------------------------------------------------------------------------------------------------


pdf('8_BarPlot_Gene_Lists_data.pdf')
gene_iter <- c( "GJA1",  "P2RY12", "PLP1",  "LY6C1","SCN7A", "C1QL1")
jitter_ics(object = Combined,target = "Cell_Type",features = gene_iter,type = "bar",color_list = color_list,organism="mouse",normalized="counts")
gene_iter <-c("CCDC153", "KCNJ13", "MEG3" ,"SNHG11", "CALD1", "VTN", 'NOTCH3')
jitter_ics(object = Combined,target = "Cell_Type",features = gene_iter,type = "bar",color_list = color_list,organism="mouse",normalized="counts")
dev.off()





features=c("PLPP3", "SLC1A2", "GJA1", "AQP4",
           "HEXB", "P2RY12", "CX3CR1", "SIGLECH",
           "PLP1", "MBP", "MOBP", "TRF",
           "LY6C1", "CLDN5", "PLTP", "PECAM1",
           "SCN7A", "MAP2", "C1QL1", "PDGFRA",
           "CCDC153", "TMEM212", "DYNLRB2", "RSPH4A",
           "KCNJ13", "FOLR1", "CLIC6", "KL",
           "MEG3", "SNHG11", "NDRG4", "ENO2",
           "CALD1", "VTN", "NOTCH3", "SNAP25")

features=c("GJA1", "P2RY12",
           "PLP1", "LY6C1","SCN7A", "C1QL1","CCDC153","KCNJ13","MEG3", "SNHG11","CALD1", "VTN", "NOTCH3")

debugonce(ics_scanpy)
ics_scanpy(object=Combined,
           group.by="Cluster",
           features=features,
           assay="RNA",
           method="heat")


write.table(cbind(Combined@meta.data,
                  Combined@reductions$umap@cell.embeddings,
                  Combined@reductions$tsne@cell.embeddings,
                  t(Combined@assays$RNA@counts)),"Dataset_GEM_Umap_Metadata.tsv",sep="\t")


write.table(cbind(Micro_Sub2@meta.data,
                  Micro_Sub2@reductions$umap@cell.embeddings,
                  Micro_Sub2@reductions$tsne@cell.embeddings,
                  t(Micro_Sub2@assays$RNA@counts)),"Microglia_GEM_Umap_Metadata.tsv",sep="\t")
