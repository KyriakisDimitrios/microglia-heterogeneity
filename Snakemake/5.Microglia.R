
source("C:/Users/dimitrios.kyriakis/Desktop/Rescue_MBSYN/Snakemake/5.Functions.R")

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
library(VGAM)
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







# =============================== Microglia Analysis Integraded  ===============================
# =============================== Microglia Analysis Integraded  ===============================
# ======================================== Microglia Analysis  =================================
dir.create("5.Microglia_Analysis")
setwd("5.Microglia_Analysis")
color_clust_microglia <- c("#42858C" ,rgb(128, 128, 192, maxColorValue = 255),
    rgb(255,128,192, maxColorValue = 255),
    #rgb(128,0,64, maxColorValue = 255),
    rgb(128,0,255, maxColorValue = 255),"blue")
color_cells_microglia <- c(rgb(128,0,128, maxColorValue = 255))
color_list_microglia <- list(condition=color_cond,Cluster=color_clust_microglia,Cell_Type=color_cells_microglia,State=color_clust_microglia)




Micro_Sub <- object_subset(object=Combined,Conditions=NULL,Clusters=NULL,Cell_Types=c("Microglia") )
DefaultAssay(Micro_Sub) <- "integrated"
Micro_Sub <- reduce_dim(Micro_Sub,project=project,resolution=c(1.2))$Combined


pdf("Umap.pdf",width=8,height=8)
DimPlot(Micro_Sub,group.by = "condition",reduction = "umap",cols=color_list_microglia$condition)
DimPlot(Micro_Sub,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)
DimPlot(Micro_Sub,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)
dev.off()

pdf("Tsne.pdf",width=8,height=8)
DimPlot(Micro_Sub,group.by = "condition",reduction = "tsne",cols=color_list_microglia$condition)
DimPlot(Micro_Sub,group.by = "Cluster",reduction = "tsne",cols=color_list_microglia$Cluster)
DimPlot(Micro_Sub,group.by = "Cluster",reduction = "tsne",cols=color_list_microglia$Cluster)
dev.off()



# =================================== DF Clusters =================================================
dir.create("DF_Clusters")
setwd("DF_Clusters")
return_fun <- df_genes(Micro_Sub,"Cluster",top_n=30,logfc.threshold=0 ,n_cores=4,latent.vars=c("nCount_RNA"),only.pos=TRUE)
setwd("../")





# =========================================================== RE-ANALYSIS ==================================================================
# =========================================================== RE-ANALYSIS ==================================================================
# =========================================================== RE-ANALYSIS ==================================================================
dir.create("Reanalysis_Counts")
setwd("Reanalysis_Counts")
#Remove Oligos from Microglia
Micro_Sub2 <- object_subset(object=Micro_Sub,Conditions=NULL,Clusters=c(1,2,3),Cell_Types=NULL)
DefaultAssay(Micro_Sub2) <- "RNA"
Micro_Sub2 <- reduce_dim(Micro_Sub2,project=project,resolution=c(1.1))$Combined
DimPlot(Micro_Sub2,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)


temp <- Micro_Sub2$Cluster  
temp[temp==1]
gene_list <- c("P2RY12","GJA1","MBP","LY6C2","TTR","ENPP2","PDGFRA","MEG3")
plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster")

pdf("Umap.pdf",width=8,height=8)
DimPlot(Micro_Sub2,group.by = "condition",reduction = "umap",cols=color_list_microglia$condition)
DimPlot(Micro_Sub2,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)
DimPlot(Micro_Sub2,group.by = "Cluster",reduction = "umap",cols=color_list_microglia$Cluster)
dev.off()

pdf("Tsne.pdf",width=8,height=8)
DimPlot(Micro_Sub2,group.by = "condition",reduction = "tsne",cols=color_list_microglia$condition)
DimPlot(Micro_Sub2,group.by = "Cluster",reduction = "tsne",cols=color_list_microglia$Cluster)
DimPlot(Micro_Sub2,group.by = "Cluster",reduction = "tsne",cols=color_list_microglia$Cluster)
dev.off()



# =================================== DF Clusters =================================================
dir.create("DF_Clusters")
setwd("DF_Clusters")
return_fun_cl_micro2 <- df_genes(Micro_Sub2,"Cluster",top_n=30,
    logfc.threshold=0 ,n_cores=4,latent.vars=c("nCount_RNA"),
    only.pos=TRUE)
color_list <- color_list_microglia
gene_list <- c("LGMN","MARCKS","P2RY12","FCRLS","CSF1R","SERINC3")
plot_bar_cells(gene_list,object=Micro_Sub2,save=TRUE,target="Cluster",title="Homeostatic")
setwd("../")



# ================================= MONOCLE =============================
# ================================= MONOCLE =============================
# ================================= MONOCLE =============================
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

pdf('DF_VLN_Monocle.pdf')
DefaultAssay(Micro_Sub2)<-"RNA"
plots <-VlnPlot(object = Micro_Sub2,features = top10_genes,group.by = "Cluster",ncol = 1,combine=FALSE)
plots <- lapply(
  X = plots,
  FUN = function(p) p + ggplot2::scale_fill_manual(values = color_list_microglia[["Cluster"]])
)
#CombinePlots(plots = plots, legend = 'right',ncol =1 )
plots
DefaultAssay(Micro_Sub2)<-"integrated"
dev.off()


setwd("../") # End Monocle











dir.create("Paper_Figures")
setwd("Paper_Figures")
color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
color_cells <- c(rgb(255,0,0, maxColorValue = 255),rgb(0,0,255, maxColorValue = 255),rgb(0,0,160, maxColorValue = 255),rgb(255,100,177, maxColorValue = 255),rgb(128,0,128, maxColorValue = 255),rgb(128,64,0, maxColorValue = 255),rgb(0,255,0, maxColorValue = 255))
color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)


p1 <-DimPlot(Combined,reduction="tsne",group.by="condition",cols=color_list$condition)
p2 <-DimPlot(Combined,reduction="tsne",group.by="Cluster",cols=color_list$Cluster)
p3 <-DimPlot(Combined,reduction="tsne",group.by="Cell_Type",cols=color_list$Cell_Type)

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



color_clust_microglia <- c(rgb(128, 128, 192, maxColorValue = 255),
    rgb(255,128,192, maxColorValue = 255),
    rgb(128,0,64, maxColorValue = 255),
    rgb(128,0,255, maxColorValue = 255))
color_cells_microglia <- c(rgb(128,0,128, maxColorValue = 255))
color_list_microglia <- list(condition=color_cond,Cluster=color_clust_microglia,Cell_Type=color_cells_microglia,State=color_clust_microglia)




dir.create("Micro_Sub2")
setwd("Micro_Sub2")

p1 <-DimPlot(Micro_Sub2,reduction="tsne",group.by="condition",cols=color_list$condition)
p2 <-DimPlot(Micro_Sub2,reduction="tsne",group.by="Cluster",cols=color_list$Cluster)
p3 <-DimPlot(Micro_Sub2,reduction="tsne",group.by="Cell_Type",cols=color_list$Cell_Type)

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

p1 <-DimPlot(Micro_Sub2,reduction="umap",group.by="condition",cols=color_list$condition)
p2 <-DimPlot(Micro_Sub2,reduction="umap",group.by="Cluster",cols=color_list$Cluster)
p3 <-DimPlot(Micro_Sub2,reduction="umap",group.by="Cell_Type",cols=color_list$Cell_Type)

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

color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
color_clust_microglia <- c("#42858C",rgb(128, 128, 192, maxColorValue = 255),
                           rgb(128,0,255, maxColorValue = 255),
                           rgb(255,128,192, maxColorValue = 255),
                           rgb(128,0,64, maxColorValue = 255))

p1 <- DimPlot(Micro_Sub2,group.by = "Cluster",cols= color_clust_microglia )


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


    
features <- c("SOCS3","CD14","GPR84","LYZ1","CASP4","FTH1","ICAM1","IL1B","ADAMTS1","CD83","CCL4","NFKBIZ")
p3<-DotPlot(object = Micro_Sub2, features = features,group.by = "micsub")+ theme(axis.text.x = element_text(angle = 45, hjust=1),axis.text.y = element_text(angle =30, hjust=1))

features <- c("H2.AA","H2.AB1","CD74")
p4 <- DotPlot(object = Micro_Sub2, features = features,group.by = "micsub")+ theme(axis.text.x = element_text(angle = 45, hjust=1),axis.text.y = element_text(angle =30, hjust=1))

features <- c("HEXB","CX3CR1","P2RY12","C1QA","FCRLS")
p5 <- DotPlot(object = Micro_Sub2, features = features,group.by = "")+ theme(axis.text.x = element_text(angle = 45, hjust=1),axis.text.y = element_text(angle =30, hjust=1))


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


pdf("DotPlots.pdf",width=12)
p3
p4
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
