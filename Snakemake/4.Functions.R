
#' data_selection
#' @author Dimitrios Kyriakis
#' @export
data_selection <- function(project){
    if(get_os()=="windows"){
        WORKDIR=paste0("C:/Users/dimitrios.kyriakis/Desktop/PhD/Projects/",project,"/")
    }else{
        WORKDIR=paste0("/home/users/dkyriakis/PhD/Projects/",project,"/")
    }
    if(project=="MBSYN"){
        cond="WT"
        if(cond=="WT"){
            list_of_files   <- c(paste(WORKDIR,"DATA/MBSYN3_S1_DGE.txt",sep=""),
                               paste(WORKDIR,"DATA/MBSYN4_S2_DGE.txt",sep=""))
            condition_names <- c("WT mice (Midbrain)",
                               "WT mice (Striatum)")

        }else if (cond=="TG"){
            list_of_files   <- c(paste(WORKDIR,"DATA/MBSYN1_S1_DGE.txt",sep=""),
                               paste(WORKDIR,"DATA/MBSYN2_S2_DGE.txt",sep=""))
            condition_names <- c("TG a-syn (Midbrain)",
                               "TG a-syn (Striatum)")
        }else{
            list_of_files   <- c(paste(WORKDIR,"DATA/MBSYN3_S1_DGE.txt",sep=""),
                               paste(WORKDIR,"DATA/MBSYN4_S2_DGE.txt",sep=""),
                               paste(WORKDIR,"DATA/MBSYN1_S1_DGE.txt",sep=""),
                               paste(WORKDIR,"DATA/MBSYN2_S2_DGE.txt",sep=""))

            condition_names <- c("WT mice (Midbrain)",
                               "WT mice (Striatum)",
                               "TG a-syn (Midbrain)",
                               "TG a-syn (Striatum)")

        }
        file <- paste0(WORKDIR,"/Gene_Lists/Cell_types2.txt")
        organism <- "mouse"
        data_10x=FALSE
    }
    return(list("WORKDIR"=WORKDIR,"list_of_files"=list_of_files,"condition_names"=condition_names,"organism"=organism,"file"=file,"data_10x"=data_10x))
}



#' Identify Operation System.
#' @author Dimitrios Kyriakis
#' @export
#' @return Returns the Operation System
#' @examples
#' get_os()
get_os <- function(){
    sysinf <- Sys.info()
    if (!is.null(sysinf)){
        os <- sysinf['sysname']
    if (os == 'Darwin')
        os <- "osx"
    } else { ## mystery machine
        os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
        os <- "osx"
    if (grepl("linux-gnu", R.version$os))
        os <- "linux"
    }
    tolower(os)
}



#' object_identifier
#'
#' @author Dimitrios Kyriakis
#' @export
object_identifier<- function(object){
    if(class(object)=="Seurat"){
        tool="seurat"
    }else{
        tool="monocle"
    }
    return(tool)
}







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
    color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
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




#' Convert Gene list to Mouse gene list
#' @author Dimitrios Kyriakis
#' @export
#' @return Gene list as Mouse gene list
#' @examples
#' simpleCap(x)
simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    s<-paste(toupper(substring(s, 1,1)), substring(s, 2),
           sep="", collapse=" ")

    if(length(grep("Loc",s) ==0)){
       s<- toupper(s)
    }
    if(length(grep("^Ay",s) ==0)){
        s<- toupper(s)
    }
    return(s)
}







#' Differential Expression Analysis.
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: object class from Monocle/Seurat.
#' @param target: target for the differential expression analysis.
#' @param top_n: Number of genes to plot for each cluster.
#' @param logfc.threshold: Threshold for Fold Change
#'
#' @return Gene expression Matrix without Outliers
#' @examples
#' df_genes(Seurat,target,top_n,name)
df_genes <- function(object,target,top_n=15,logfc.threshold=0,n_cores=4,latent.vars=NULL, min.pct = 0,only.pos = FALSE,title="DF_Genes",p_val_adj_thres=0.05){
    tool <- object_identifier(object)
    if(tolower(tool)=="monocle"){
        all_diff_test_res_condition <- differentialGeneTest(object, fullModelFormulaStr=paste("~",target,sep=""),cores=n_cores, reducedModelFormulaStr ="~num_genes_expressed")
        # == Reduce Matrix (remove not significant genes)
        signif_genes <-all_diff_test_res_condition
        monocle_helper <- df_genes_logfc(object,target,signif_genes=signif_genes,top_n=top_n,qval_thres=0.05,fc_thres=logfc.threshold,each_cl=TRUE,n_cores=4)
        pbmc.markers <- monocle_helper$df_pval_genes
        colnames(pbmc.markers)[colnames(pbmc.markers)=="gene_short_name"] <- "gene"
        top10_genes <- monocle_helper$sig_genes

    }else{
        DefaultAssay(object) <- "RNA"
        object <- ScaleData(object)
        Idents(object) <- object[[target]]
        # find markers for every cluster compared to all remaining cells, report only the positive ones
        pbmc.markers <- FindAllMarkers(object = object,
                                       assay ="RNA",
                                       logfc.threshold=logfc.threshold,
                                       only.pos = only.pos,
                                       min.pct=min.pct,
                                       test.use = "MAST",latent.vars = latent.vars)
        qvalue <- p.adjust(pbmc.markers$p_val, method = "BH",n=dim(object@assays$RNA@counts)[1])
        pbmc.markers$qvalue <- qvalue
        top <- pbmc.markers[pbmc.markers$p_val_adj<p_val_adj_thres,]
        top10 <- top %>% group_by(cluster) %>% top_n(n = top_n, wt = abs(avg_logFC))
        top10_genes<- unique(top10$gene)
    }


    #debugonce(annotated_heat)
    annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name=title,
                  title="DF_genes",
                  ordering=target)
    if(tolower(tool)!="monocle"){
        DefaultAssay(object) <- "integrated"
    }

    write.table(pbmc.markers, file = paste0(Sys.Date(),"_TO_EXP_each_",target,"_",title,".tsv"),row.names=FALSE, na="", sep="\t")
    return(list("degs"=pbmc.markers,"top_markers"=top10_genes))
}




#' plot_bar_cells
#' @author Dimitrios Kyriakis
#' @export
plot_bar_cells <- function(gene_list,object,save=TRUE,target="Cell_Type",title=""){
  # Plot the log counts of genes from a matrix
  # object <- Cell Data
  tool <- object_identifier(object)
  plist <- lapply(gene_list,barplot_gene,object=object,target=target)
  if(save){
    figure<- ggarrange(plotlist=plist, nrow=length(gene_list), common.legend = TRUE, legend="bottom")
    annotate_figure(figure,
                    # bottom = text_grob("Cells"),
                    left = text_grob("Counts", rot = 90),
    )%>%
      ggexport(filename = paste(Sys.Date(),title,"Barplot_marker_genes.pdf",sep="_"))
    dev.off()
  }else{
    figure<- ggarrange(plotlist=plist, nrow=length(gene_list), common.legend = TRUE, legend="bottom")
    annotate_figure(figure,
                    # bottom = text_grob("Cells"),
                    left = text_grob("Counts", rot = 90),
    )
  }
}



#' barplot_gene
#' @author Dimitrios Kyriakis
#' @export
barplot_gene <- function(x,object,target="Cell_Type"){
  gene_list=x
  tool <- object_identifier(object)
  if(tolower(tool) =="seurat"){
    GE_matrix <- as.matrix(object@assays$RNA@counts)
    # range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    #
    # GE_matrix <- range01(object@assays$integrated@data)
    # GE_matrix <- object@assays$integrated@data
    # GE_matrix[GE_matrix<0] <- 0
  }else{
    GE_matrix<- data.matrix(object, rownames.force = NA)
    GE_matrix <- as.matrix(GE_matrix)
    rownames(GE_matrix) <- fData(object)$gene_short_name

  }
  if (target=="Cluster"){
    ordering<- order(object$Cluster)
    Cell_Types <- object$Cluster[ordering]
  }else{
    ordering<- order(object$Cell_Type)
    Cell_Types <- object$Cell_Type[ordering]
  }

  gene_indx <- match(x,rownames(GE_matrix))


  index_cells <- match(levels(as.factor(Cell_Types)),Cell_Types)
  # number_of_transcripts <- rowSums(GE_matrix)
  # names(number_of_transcripts) <- fData(object)$gene_short_name

  df<- data.frame(names = c(1:dim(object)[2]),y=GE_matrix[gene_indx,ordering] )

  if(x==gene_list[length(gene_list)]){
    p1 <- ggplot(data=df, aes(x=names, y=y,fill=Cell_Types)) +
      geom_bar(stat="identity")+
      ggtitle(simpleCap(tolower(x)))+xlab("")+ylab("")+
      scale_fill_manual(values = color_list[["Cell_Type"]], name = "Cell Types",na.value = "gray")+

      # TO CHnage
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    # scale_x_continuous(breaks=index_cells,labels=levels(as.factor(Cell_Typess)))+
    # scale_x_continuous(breaks=index_cells,labels=c("Astr","BlC","End","Epe","Micr","Neur","Olig"))+

    # theme(axis.text.x = element_text(angle = 50))

  }else{
    p1 <- ggplot(data=df, aes(x=names, y=y,fill=Cell_Types)) +
      geom_bar(stat="identity")+
      scale_fill_manual(values = color_list[["Cell_Type"]], name = "Cell Types",na.value = "gray")+

      ggtitle(simpleCap(tolower(x)))+xlab("")+ylab("")+

      # ggtitle("$\textit{x}$")+xlab("")+ylab("")+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
  }
  return(p1)
}


