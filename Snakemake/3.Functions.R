
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






#' cell_type_assignment.
#' @export
#' @param object: object class from Monocle/Seurat.
#' @param tab_name: The name of the Tab in the object to store the assignment
#' @param file: .file path with the markers
#'
#' @return Gene expression Matrix without Outliers
#' @examples
#' cell_type_assignment(object,tab_name,file)
cell_type_assignment <- function(object,tab_name,file,assign=TRUE){
    tool <- object_identifier(object)
    if(tolower(tool)=="monocle"){
        GE_matrix <- data.matrix(object, rownames.force = NA)
        Normilize_Vector <- pData(object)[, 'Size_Factor']
        count_matrix <- as.matrix(sweep(GE_matrix,2,Normilize_Vector,"/"))
    }else{
        count_matrix <- as.matrix(object@assays$RNA@data)
    }
    astro <-as.matrix(read.delim2(file,header=F))
    r_annot <- astro[,1]
    gene_list <- toupper(astro[,2])
    cell_type_list <- list()

    df_cell_type = data.frame(matrix(vector(), dim(object)[2]))

    for (cell_type in unique(r_annot)){
        print(cell_type)
        cell_type_list<- list(c(gene_list[grep(cell_type,r_annot)]))
        indexes <- cell_type_list[[1]]%in% rownames(object)
        gene_list_s <- cell_type_list[[1]][indexes]
        ctls <- length(cell_type_list[[1]])

        if(length(gene_list_s)==0){
            cat(yellow("Warning :")  %+% "This Cell Type skipped due to no genes found\n")
        }else if(length(gene_list_s)==1){
            df_cell_type[[cell_type]] <- as.vector((count_matrix[gene_list_s,]))
        }else{
            df_cell_type[[cell_type]] <- as.vector(colSums(count_matrix[gene_list_s,]))/length(gene_list_s)
        }

        if(length(gene_list_s)!=0){
            a <- 1:length(gene_list_s)
            n <- length(a)
            k <- 4
            splited_l <- split(a, rep(1:ceiling(n/k), each=k)[1:n])

            print_fun<-function(x) {
                if(tolower(tool)=="monocle"){
                    plot(plot_cell_clusters(object, x = 1, y = 2, markers = gene_list_s[x]))
                }else{
                    print(FeaturePlot(object = object, features = gene_list_s[x],reduction = "tsne",cols=brewer.pal(100, "OrRd"),order=T))
                }

            }

            pdf(paste(Sys.Date(),"Scat",project,cell_type,"tsne.pdf",sep="_"))
            lapply(splited_l,print_fun)
            dev.off()
        }

    }

    cell_type_vector<-as.vector(object$Cluster)
    df_cl_cell_type = data.frame(matrix(vector(),nrow = dim(df_cell_type)[2]))
    rownames(df_cl_cell_type)<- colnames(df_cell_type)
    max_num_cell_type <- list()
    for (cl in levels(object$Cluster)){
        Cluster_score <- colSums(df_cell_type[object$Cluster==cl,])
        cl_cell_type <- colnames(df_cell_type)[which.max(Cluster_score)]

        max_num <- Cluster_score[which.max(Cluster_score)]
        name_cell <- names(max_num)
        if(name_cell%in%names(max_num_cell_type)){
          if(max_num>=max_num_cell_type[name_cell]){
            max_num_cell_type[name_cell] <- max_num
            cell_type_vector[object$Cluster==cl] <- cl_cell_type
          }else{
            cell_type_vector[object$Cluster==cl] <- "Hybrid"
          }
        }else{
          max_num_cell_type[name_cell] <- max_num
          cell_type_vector[object$Cluster==cl] <- cl_cell_type
        }

        df_cl_cell_type[,cl] <- c(Cluster_score)
    }

    if (assign){
        object[[tab_name]] <- cell_type_vector


        annotated_heat(object=object,
                     row_annotation=as.vector(r_annot),
                     gene_list=gene_list,
                     gene_list_name="Markers",
                     title="Cell Type Markers",
                     ordering="Cluster")
    }

    return(list("object"=object,"r_annot"=r_annot))
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

    color_cond  <- c(rgb(0,128,0, maxColorValue = 255),rgb(255, 128,0, maxColorValue = 255))
    color_clust <- c(brewer.pal(9,"Set1")[-6],"goldenrod4","darkblue","seagreen4")
    color_cells <- c(rgb(255,0,0, maxColorValue = 255),rgb(0,0,255, maxColorValue = 255),rgb(0,0,160, maxColorValue = 255),rgb(255,100,177, maxColorValue = 255),rgb(128,0,128, maxColorValue = 255),rgb(128,64,0, maxColorValue = 255),rgb(0,255,0, maxColorValue = 255))
    color_list <- list(condition=color_cond,Cluster=color_clust,Cell_Type=color_cells,State=color_clust)

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




