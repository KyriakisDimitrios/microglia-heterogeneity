#' barplot_gene
#' @author Dimitrios Kyriakis
#' @export
barplot_gene <- function(x,object,target="Cell_Type"){
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


#' pca_elbow
#' @author Dimitrios Kyriakis
#' @export
pie_per_CellType<-function(object,target) {
  nums <- table(object[[target]])
  pie_data<- data.frame(value = as.vector(nums),
                        Groups = paste0(names(nums)," ",round(as.vector(nums)/sum(as.vector(nums)) *100 ,1),"%")   ) %>%
    mutate(Group=factor(Groups, levels=Groups),
           cumulative = cumsum(value),
           midpoint = cumulative - value /2 ,
           label = paste0(Groups," ",round(value/sum(value) *100 ,1),"%"),
           label2 = paste0(round(value/sum(value) *100),"%"))
  pie<-ggplot(pie_data,
              aes(x= factor(1),weight=value, fill = Groups)) +
    geom_bar(width=1,position="stack") +
    coord_polar(theta="y")+
    scale_color_manual(values = color_list["Clusters"])

  pdf(paste(Sys.Date(),"_Pie_CellType_",project,".pdf",sep=""))
  print(pie)
  dev.off()
}



#' pca_elbow
#' @author Dimitrios Kyriakis
#' @export
pca_elbow <- function(object,project){
    tool <- object_identifier(object)
    pdf(paste(Sys.Date(),"PCA-Elbow_",project,".pdf"))
    if(tool=="monocle"){
        plot(plot_pc_variance_explained(object, return_all = F))
    }else{
        plot(ElbowPlot(object = object))
    }
    dev.off()
}








#' Plot scatter plot of the dataset in 2 dimension tSNE projection
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: object class from Monocle/Seurat.
#' @param target: the condition to color by (Default= condition)
#' @param reduction: Reduction method
#' @param only_cond: If you want to color only specific condition/cluster
#' @param save: Save as pdf or plot
#' @param leg_pos: Legend position (Default= bottom)
#' @param extra_title: Extra title
#' @param ncol: number of columns in legend
#'
#' @return Plot scatter plot
plot_cells <- function(object,target="condition",leg_pos="bottom",save=TRUE,extra_title=NA,ncol=1,only_cond=NA,reduction="tsne",p3D=FALSE){
    tool <- object_identifier(object)
    colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
    color_cond  <- c(brewer.pal(5,"Dark2"),"black","gray","magenta4","seagreen4")
    color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
    color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
    color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)
    #color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

    color_list_temp <- color_list

    if (tolower(tool)=="monocle"){
        # ========== Option to plot the pca
        label_xy <- reduction
        if(toupper(reduction)=="PCA"){
            object@reducedDimA <- object@auxClusteringData[["tSNE"]]$reduced_dimension[1:3,]
            Vars <- object@auxClusteringData[["tSNE"]]$variance_explained[1:3]
            label_x <- paste("PC1", " (",round(Vars[1]*100,2),"% explained var.)",sep="")
            label_y <- paste("PC2", " (",round(Vars[2]*100,2),"% explained var.)",sep="")
        }else if(reduction=="dm"){
            object@reducedDimA <- object@reducedDimS
            label_x="Component_1"
            label_y="Component_2"
        }else{
            p3D=FALSE
            label_x="t-SNE1"
            label_y="t-SNE2"
        }
        # -----------------------------

        if (!is.na(only_cond)){
            color_list_temp[[target]] <-  color_list_temp[[target]][levels(object[[target]])%in%only_cond]
            object[[target]][!object[[target]] %in% only_cond] <- NA
        }
        object[[target]] <- factor(object[[target]])
        title<-paste(paste("cond-",length(levels(object$condition)),sep=""),paste("clust-",length(levels(object$Cluster)),sep=""),sep="_")

        # Figure Title
        fig_title <- paste(Sys.Date(),"_",dataset,"_Scat_",target,"_",title,sep="")
        if(!is.na(extra_title)){fig_title<-paste(fig_title,extra_title,sep="_")}
        if(!is.na(only_cond)){fig_title<-paste(fig_title,only_cond,sep="_")}

        if (save){ pdf(paste(fig_title,".pdf",sep=""))}

        if(p3D){
            data_tot <- data.frame(PC1=object@reducedDimA[1,],PC2 = object@reducedDimA[2,],PC3=object@reducedDimA[3,],counts = object$Total_mRNA)
            label_z <- paste("PC3", " (",round(Vars[3]*100,2),"% explained var.)",sep="")

            p1 <- plot_ly(data_tot, x = ~PC1, y = ~PC2, z = ~PC3,size=~log10(counts),color = ~object[[target]], showscale = TRUE,colorscale=color_list_temp[[target]])%>%
              layout(title = 'Quality Control of Clusters based on Total mRNA',
                     scene = list(xaxis = list(title = label_x),
                                  yaxis = list(title = label_y),
                                  zaxis = list(title = label_z)))


        }else{
            p1<- plot_cell_clusters(object,1,2,color_by =target)+theme(legend.position=leg_pos)+
          guides(col=guide_legend(ncol=ncol))+
          scale_color_manual(values = color_list_temp[[target]], name = target,na.value = "gray")+
          xlab(label_x)+ylab(label_y)
        }
    }else{
        p1 <- DimPlot(object = object, reduction = reduction, group.by = target,cols = color_list_temp[[target]])
    }



    if(!is.na(only_cond)){print(p1+theme(legend.position = "none")+ggtitle(only_cond))} else{print(p1)}
    if (save){dev.off()}
    return(p1)
}
# ------------------------------------------------------- #




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
        object_rownames <-rownames(object)
        GE_matrix<- GetAssayData(object = object, slot = "scale.data")
        GE_matrix <- as.matrix(GE_matrix)
    }
    disp.min <- -2.5
    disp.max <- 2.5
    test0 <- MinMax(data = GE_matrix, min = disp.min, max = disp.max)

    list_ind1 <- match(gene_list,object_rownames)

    na_index <- is.na(list_ind1)
    r_annotation<-data.frame(gene_list_name=r_annotation[!na_index])
    gene_list <- gene_list[!na_index]
    list_ind1<- list_ind1[!na_index]

    # Subset Top genes
    Top_DF_Cond <- test0[list_ind1,]
    rownames(Top_DF_Cond)<-gene_list


    colormap_d<- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928',"black","gray")
    color_cond  <- c(brewer.pal(8,"Dark2"),"black","gray","magenta4","seagreen4")
    color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
    color_cells <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
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
                # color=CustomPalette(low = "#000000", high = "#ff0000", k = 100),
                cexRow=3.5,
                color=CustomPalette(low = "#ff0905", mid = "#000000" , high = "#ffff05", k = 100),
                #CustomPalette(low = "#000000" , high = "#db1200", mid = "#400459", k = 50),
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
                color=CustomPalette(low = "#ff0905", mid = "#000000" , high = "#ffff05", k = 100),#rev(brewer.pal(name = "RdYlBu",n =11)),#,jcolors(palette = "pal12"),
                labRow = NULL,
                main=paste("Gene marker for ",title),
                filename=paste(Sys.Date(),"_Heat_",title,"_target-",ordering,".pdf",sep=""),width = 13, height = 10)
    }

}

#' plot_tot_mRNA
#' @author Dimitrios Kyriakis
#' @export
plot_tot_mRNA <- function(object,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE){
    tool <- object_identifier(object)
    if(tolower(tool)=="monocle"){
        # Plot the log counts of total counts ~= Total mRNA
        # object <- Cell Data

        # ========== Option to plot the pca
        label_xy <- reduce

        if(toupper(reduce)=="PCA"){
            object@reducedDimA <- object@auxClusteringData[["tSNE"]]$reduced_dimension[1:3,]
            Vars <- object@auxClusteringData[["tSNE"]]$variance_explained[1:3]
            label_x <- paste("PC1", " (",round(Vars[1]*100,2),"% explained var.)",sep="")
            label_y <- paste("PC2", " (",round(Vars[2]*100,2),"% explained var.)",sep="")
        }else{
            p3D=FALSE
            label_x="t-SNE1"
            label_y="t-SNE2"
        }
        # -----------------------------
        data_tot <- data.frame(tSNE1=object@reducedDimA[1,],tSNE2 = object@reducedDimA[2,],counts = object$nCount_RNA)

        if(p3D){
            data_tot <- data.frame(PC1=object@reducedDimA[1,],PC2 = object@reducedDimA[2,],PC3=object@reducedDimA[3,],counts = object$nCount_RNA)
            label_z <- paste("PC3", " (",round(Vars[3]*100,2),"% explained var.)",sep="")


            p1<-plot_ly(data_tot, x = ~PC1, y = ~PC2, z = ~PC3,size=~log1p(counts),color = ~log1p(counts), showscale = TRUE,colorscale=jcolors("pal12"))%>%
              layout(title = 'Quality Control of Clusters based on Total mRNA',
                     scene = list(xaxis = list(title = label_x),
                                  yaxis = list(title = label_y),
                                  zaxis = list(title = label_z)))

        }else{
            p1 <- ggplot(data_tot, aes(x = tSNE1, y = tSNE2,color=log1p(counts))) +
              geom_point(mapping = aes(size = counts^2)) +
              scale_color_jcolors_contin("pal12")+
              xlab(label_x)+ylab(label_y)+
              ggtitle("Quality Control of Clusters based on Total mRNA")+
              scale_size(name   = "Total mRNA",
                         breaks = fivenum(data_tot$counts)^2,
                         labels = fivenum(data_tot$counts))
        }


    }else{
        p1 <- FeaturePlot(object, features = c("nCount_RNA"),
            reduction = "tsne",
            pt.size=1.2,
            cols=jcolors(palette = "pal12"))
    }
    if(save){
        filename <- paste(Sys.Date(),"_Scat_",title,"QC_Total_mRNA.pdf",sep="")
        pdf(filename)
        print(p1)
        dev.off()
    }else{
        print(p1)
    }
    if(tiff){
        filename <- paste(Sys.Date(),"_Scat_",title,"QC_Total_mRNA.tif",sep="")
        ggsave(filename,p1,"tiff")
    }
    return(p1)
}

#' plot_nFeatures
#' @author Dimitrios Kyriakis
#' @export
plot_nFeatures <- function(object,title="",save=TRUE,tiff=FALSE,reduce="t-SNE",p3D=FALSE){
    tool <- object_identifier(object)
    if(tolower(tool)=="monocle"){
        # Plot the log counts of nFeatures
        # object <- Cell Data
        # ========== Option to plot the pca
        label_xy <- reduce
        if(toupper(reduce)=="PCA"){
            object@reducedDimA <- object@auxClusteringData[["tSNE"]]$reduced_dimension[1:3,]
            Vars <- object@auxClusteringData[["tSNE"]]$variance_explained[1:3]
            label_x <- paste("PC1", " (",round(Vars[1]*100,2),"% explained var.)",sep="")
            label_y <- paste("PC2", " (",round(Vars[2]*100,2),"% explained var.)",sep="")

        }else{
            p3D=FALSE
            label_x="t-SNE1"
            label_y="t-SNE2"
        }
        # -----------------------------

        data_tot <- data.frame(tSNE1=object@reducedDimA[1,],tSNE2 = object@reducedDimA[2,],counts = object$num_genes_expressed)
        if(p3D){
            data_tot <- data.frame(PC1=object@reducedDimA[1,],PC2 = object@reducedDimA[2,],PC3=object@reducedDimA[3,],counts = object$num_genes_expressed)
            label_z <- paste("PC3", " (",round(Vars[3]*100,2),"% explained var.)",sep="")

            p1<-plot_ly(data_tot, x = ~PC1, y = ~PC2, z = ~PC3,size=~log1p(counts),color = ~log1p(counts), showscale = TRUE,colorscale=jcolors("pal12"))%>%
              layout(title = 'Quality Control of Clusters based on nFeatures',
                     scene = list(xaxis = list(title = label_x),
                                  yaxis = list(title = label_y),
                                  zaxis = list(title = label_z)))

        }else{

        p1 <- ggplot(data_tot, aes(x = tSNE1, y = tSNE2,color=log1p(counts))) +
          geom_point(mapping = aes(size = counts^2)) +
          scale_color_jcolors_contin("pal12")+
          xlab(label_x)+ylab(label_y)+
          ggtitle("Quality Control of Clusters based on nFeatures")+
          scale_size(name   = "nFeatures",
                     breaks = fivenum(data_tot$counts)^2,
                     labels = fivenum(data_tot$counts))
        }

    }else{
        p1 <- FeaturePlot(object, features = c("nFeature_RNA"),
            reduction = "tsne",
            pt.size=1.2,
            cols=jcolors(palette = "pal12"))
    }
    if(save){
        filename <- paste(Sys.Date(),"_Scat_",title,"QC_nFeatures.pdf",sep="")
        pdf(filename)
        print(p1)
        dev.off()
        #ggsave(filename,p1,"pdf")
    }else{
        print(p1)
    }
    if(tiff){
        filename <- paste(Sys.Date(),"_Scat_",title,"QC_nFeatures.tif",sep="")
        ggsave(filename,p1,"tiff")
    }
    return(p1)
}
