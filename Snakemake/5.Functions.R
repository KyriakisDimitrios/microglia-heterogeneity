
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
    color_clust_microglia <- c("#42858C" ,rgb(128, 128, 192, maxColorValue = 255),
    rgb(255,128,192, maxColorValue = 255),
    #rgb(128,0,64, maxColorValue = 255),
    rgb(128,0,255, maxColorValue = 255),"blue")
    color_cells_microglia <- c(rgb(128,0,128, maxColorValue = 255))
    color_list_microglia <- list(condition=color_cond,Cluster=color_clust_microglia,Cell_Type=color_clust_microglia,State=color_clust_microglia)
    color_list <-color_list_microglia 

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






#' Subsetting  Object
#' @author Dimitrios Kyriakis
#'
#' @param object : CellDataSet class from Monocle.
#' @param Conditions : Vector with the Conditions to keep.
#' @param Clusters : Vector with the Clusters to keep.
#' @param Cell_Types : Vector with the Cell_Types to keep.
#' @param tool : Monocle/Seurat
#' @export
#'
#' @return A subset of the Monocle Object.
#' @examples
#' object_subset(object,Conditions,Clusters,tool="monocle")
object_subset <- function(object,Conditions=NULL,Clusters=NULL,Cell_Types=NULL,re_project=TRUE){
    tool <- object_identifier(object)
    if(is.null(Conditions)){
        Conditions <- as.vector(unique(object$condition))
    }
    if(is.null(Clusters)){
        Clusters <- as.vector(unique(object$Cluster))
    }
    if(is.null(Cell_Types)){
        Cell_Types <- as.vector(unique(object$Cell_Type))
    }

    if(tolower(tool)=="monocle"){
        indexes <- object$Cluster%in%Clusters & object$condition%in%Conditions & object$Cell_Type%in%Cell_Types
        sub_CellType <- object[,indexes]
        print(dim(sub_CellType))
        sub_CellType@reducedDimA <- sub_CellType@reducedDimA[,indexes]
        # Re-factor the conditions
        sub_CellType$condition <- factor(sub_CellType$condition)
        sub_CellType$Cluster <- factor(sub_CellType$Cluster)
        sub_CellType$Cell_Type <- factor(sub_CellType$Cell_Type)
    }else{
        # if(!is.null(Conditions)){
        #     Idents(object) <- object[["condition"]]
        #     sub_CellType <- SubsetData(object = object, ident.use = Conditions)
        # }else{
        #     sub_CellType <- object
        # }
        # if(!is.null(Clusters)){
        #     Idents(sub_CellType) <- sub_CellType[["Cluster"]]
        #     sub_CellType <- SubsetData(object = sub_CellType, ident.use = Clusters)
        # }
        # if(!is.null(Cell_Types)){
            # Idents(object) <- object[["Cell_Type"]]
            # sub_CellType <- SubsetData(object = object, ident.use = Cell_Types)
        # }
        cells_index <- object$Cell_Type%in%Cell_Types & object$condition%in%Conditions &  object$Cluster%in%Clusters
        sub_CellType<- subset(x=object,cells=colnames(object)[cells_index])
    }
    # ============== Re_Project ==================
    if(re_project){
      sub_CellType <- reduce_dim(sub_CellType,project=project)$Combined
    }
    # --------------------------------------------
  return(sub_CellType)
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






#' Calculates how many Principal Components contain N\% of the variation of the data.
#'
#' @author Dimitrios Kyriakis
#' @export
#' @param CDSC: CellDataSet class from Monocle.
#' @param threshold: Expain variance.
#' @return The number of PCs expaining the N\% of the variation of the data.
#' @examples calc_num_pc <- function(CDSC,0.95).
calc_num_pc <- function(object,threshold){
  tool <- object_identifier(object)
  if(tolower(tool)=="monocle"){
    pc<-plot_pc_variance_explained(object, return_all = F)
    listaa <-pc$plot_env$irlba_res$sdev/100
  }else{
    x = Stdev(object = object, reduction = "pca")
    listaa<-x/100
  }
  z<-0
  for (i  in 1:length(listaa)) {
    z<- z +listaa[i]
    if (z>threshold){
      break}
  }
  cat(cyan("Explained Variance : ")  %+% paste(threshold, "% of the variance is explained by ", i ," Principal Components\n",sep=""))
  return(i)
}








#' reduce_dim
#' @author Dimitrios Kyriakis
#' @export
reduce_dim<-function(object,project,resolution=FALSE){
    tool <- object_identifier(object)
    if(tolower(tool)=="monocle"){

        object <- detectGenes(object, min_expr = 0.1)
        # Filter genes based on average expression level,
        # and we can additionally select genes that are unusually variable across cells.
        # These genes tend to be highly informative about cell state.
        # Dispersion analysis
        disp_table <- dispersionTable(object)
        unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)

        # The setOrderingFilter function marks genes that will be used for clustering in subsequent calls to clusterCells,
        # although we will be able to provide other lists of genes if we want.
        object <- setOrderingFilter(object, unsup_clustering_genes$gene_id)

    }else{
        # object <- NormalizeData(object)
        object <- FindVariableFeatures(object = object, selection.method = "vst", nfeatures = 5000)
        all.genes <- rownames(object)
        object <- ScaleData(object, features = all.genes)
        #object <- ScaleData(object = object, verbose = FALSE)
        # object <- SCTransform(object, vars.to.regress = "percent.mito", verbose = FALSE)

        object <- RunPCA(object = object,verbose = FALSE)
    }

    pca_elbow(object,project)
    num_dim <- calc_num_pc(object=object,0.95)

    if(tolower(tool)=="monocle"){
        object <- reduceDimension(object,
          max_components=2,
          num_dim = num_dim,
          reduction_method = 'tSNE',
          residualModelFormulaStr="~num_genes_expressed", verbose = T, check_duplicates=TRUE)

    }else{
        # t-SNE and Clustering
        object <- RunUMAP(object = object, reduction = "pca", dims = 1:num_dim)
        object <- RunTSNE(object = object, reduction = "pca", dims = 1:num_dim)
        object <- FindNeighbors(object = object, reduction = "pca", dims = 1:num_dim)
    }

    # =================== Optimal Clusters ======================================
    optimal_output<- optimal_clusters(object,k.max=10,save=TRUE,resolution=resolution)
    object <- optimal_output$object
    opt_num <- optimal_output$opt_num
    sil_scor <- optimal_output$sil_scor
    # ---------------------------------------------------------------------------

    return(list("Combined"=object,"Num_dim"=num_dim))
}



#' Calculate the optimal number of Clusters.
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: CellDataSet class from Monocle.
#' @param k.max: Number of Clusters to check.
#' @param save: Save image.
#' @return CellDataSet class from Monocle.
#' @return Optimal number of CLusters.
#' @examples
#' monocle_optimal_clust(object,k.max=10,save=TRUE)
optimal_clusters<- function(object,k.max=10,save=TRUE,resolution=FALSE){

    tool <- object_identifier(object)
    mean_score<-c()
    Clusts_num<-c()
    sd_score<-c()
    if(tolower(tool)=="monocle"){
        dist_mat <- dist(t(as.matrix(object@reducedDimA)))
        for (i in  c(2:k.max)){
            object <-clusterCells(object,num_clusters=i)
            if(length(levels(object$Cluster))<i){
                object <- clusterCells(object,num_clusters=i+1, verbose = F)
            }
            sil_scor<-silhouette(as.numeric(object$Cluster),dist_mat)
            mean_score[i-1]=summary(sil_scor)$si.summary[4]
            wt_size <- summary(sil_scor)$clus.sizes/sum(summary(sil_scor)$clus.sizes)
            # xm <- weighted.mean(summary(sil_scor)$clus.avg.widths, wt_size)
            xm<-mean_score[i-1]
            v <- sum(wt_size * (summary(sil_scor)$clus.avg.widths - xm)^2)/sum(wt_size)
            sd_score[i-1]=sqrt(v)
            Clusts_num[i-1]=length(levels(object$Cluster))
        }

        silhouette_df<- data.frame(Clusters=Clusts_num,Mean_Silhouette=mean_score)
        num_clusters<- which.max(mean_score)+1

        p1 <- ggplot(data=silhouette_df,aes(x=Clusters,y=Mean_Silhouette))+
            geom_point()+ geom_line()+
            # geom_pointrange(aes(ymin=Mean_Silhouette-sd_score, ymax=Mean_Silhouette+sd_score))+
            scale_x_continuous(breaks=seq(2, 15, 1))+
            geom_vline(xintercept = num_clusters, linetype="dashed", color = "red")+
            ggtitle("Optimal number of clusters")+ylab("Mean Silhouette score")+xlab("Number of Clusters")

        cat(cyan("Clustering : ")  %+% paste("The optimal number of clusters based on Mean Silhouette score is ",num_clusters,".\n"))


        object <- clusterCells(object,num_clusters=num_clusters)
        if(length(levels(object$Cluster))!= num_clusters){object <- clusterCells(object,num_clusters=num_clusters+1)}
        sil_scor<-silhouette(as.numeric(object$Cluster),dist_mat)
        optimal <- num_clusters

    }else{
        if (!resolution) {
            resolution <- c(0.5,0.6,0.7,0.8,0.9,1.0)
        }

        for (i in  c(1:length(resolution))){
            res <- resolution[i]
            object <- FindClusters(object, resolution = res)
            object$Cluster <- Idents(object = object)
            projection_mat <- Embeddings(object = object$pca)

            # Calc Silhouette
            sil_scor<-silhouette(as.numeric(object$Cluster),dist(projection_mat))
            mean_score[i]=summary(sil_scor)$si.summary[4]
            wt_size <- summary(sil_scor)$clus.sizes/sum(summary(sil_scor)$clus.sizes)
            # xm <- weighted.mean(summary(sil_scor)$clus.avg.widths, wt_size)
            xm<-mean_score[i]
            v <- sum(wt_size * (summary(sil_scor)$clus.avg.widths - xm)^2)/sum(wt_size)
            sd_score[i]=sqrt(v)
            Clusts_num[i]=length(levels(object$Cluster))

        }

        silhouette_df<- data.frame(resolution=resolution,Mean_Silhouette=mean_score)
        opt_res<- resolution[which.max(mean_score)]

        p1 <- ggplot(data=silhouette_df,aes(x=resolution,y=Mean_Silhouette))+
            geom_point()+ geom_line()+
            # geom_pointrange(aes(ymin=Mean_Silhouette-sd_score, ymax=Mean_Silhouette+sd_score))+
            scale_x_continuous()+
            geom_vline(xintercept = opt_res, linetype="dashed", color = "red")+
            ggtitle("Optimal number of clusters")+ylab("Mean Silhouette score")+xlab("Resolution")

        object <- FindClusters(object, resolution = opt_res)
        # object$Cluster <- Idents(object = object)
        object$Cluster <- as.factor(as.numeric(object$seurat_clusters))
        names(object$Cluster) <- colnames(object)

        projection_mat <- Embeddings(object = object$umap)
        sil_scor<-silhouette(as.numeric(object$Cluster),dist(projection_mat))

        print_text <- paste("The optimal Resolution based on Mean Silhouette score is ",opt_res)
        print(paste(rep("#",nchar(print_text)),collapse=""))
        print(print_text)
        print(paste(rep("-",nchar(print_text)),collapse=""))
        optimal <- opt_res
    }

    # ====================== BARPLOT CELL BELONGS TO WHICH CLUSTERS ====================================
    data <- data.frame(cbind(factor(object$Cluster),paste(" ",object$condition)))
    colnames(data)<- c("Cluster","Condition")
    pdf(paste(Sys.Date(),'Barplot_num-Cond_per_Cluster.pdf',sep="_"))
    print(ggplot(data, aes(Cluster)) +  geom_bar(aes(fill = Condition))+
              guides(col=guide_legend(ncol=2,))+
              scale_fill_brewer(palette="Dark2")+
              theme(legend.position="bottom") )
    dev.off()
    # -------------------------------------------------------------------------------------------------

    num_clusters <- length(unique(object$Cluster))

    if(save){
        pdf(paste(Sys.Date(),"pd_optimal_cluster.pdf",sep="_"))
        print(p1)
        dev.off()
        pdf(paste(Sys.Date(),"pd_silhouette.pdf",sep="_"))
        plot(sil_scor,col = color_list$Cluster[1:num_clusters])
        dev.off()
    } else{ print(p1)}

    return(list(object=object,opt_num=optimal,sil_scor=sil_scor))
}









#' seurat_to_monocle
#'
#' @author Dimitrios Kyriakis
#' @export
seurat_to_monocle<- function(object){
    mat <- as.matrix(object@assays$RNA@counts)

    data <- as(mat,"sparseMatrix")
    pd<-new("AnnotatedDataFrame",data=object@meta.data)
    pd[["num_genes_expressed"]]<-as.vector(object$nFeature_RNA)
    pd[["Total_mRNAs"]]<-as.vector(object$nCount_RNA)

    pca <- as.matrix(object@reductions$pca@cell.embeddings)
    tsne <- as.matrix(object@reductions$tsne@cell.embeddings)
    umap <- as.matrix(object@reductions$umap@cell.embeddings)

    fData <- data.frame(gene_short_name=row.names(data),row.names=row.names(data))
    fd <- new("AnnotatedDataFrame",data=fData)
    # Construct Monocle Object
    object <- newCellDataSet(data,phenoData=pd,featureData=fd,lowerDetectionLimit=0.5,
            expressionFamily=negbinomial.size())
    object <- estimateSizeFactors(object)
    #this then for dipersion!
    object <- estimateDispersions(object)
    object@reducedDimA<-t(pca)
    return(object)
}



#' Differential expression analysis between Conditions and Clusters.
#' @author Dimitrios Kyriakis
#' @export
#'
#' @family Differential Expression Analysis
#' @param object : CellDataSet class from Monocle.
#' @param target : Target for the differential expression Cluster/Condition/State.
#' @param colorby : Colored by Cluster/Condition/State.
#' @param num_thres : Number of genes to plot in Heatmap. If each_cl is True is num of genes per cluster.
#' @param normalize : Normalization of the matrix "Monocle"/"Total_mRNA" (Default Monocle).
#' @param scale : Scale the values zscore/minmax (Default=zscore).
#' @param alter : Alternative Color to the heatmap.
#' @param qval_thres : Q-value Threshold.
#' @param fc_thres : FoldChange threshold.
#' @param each_cl : Top for each cluster.
#' @param n_cores : Number of cores (Cores=1 for Windows OS).
#' @param batch_formula : Residual Formula (Default "~num_genes_expressed").

#' @return A Hist/Density plot with the average score of the gene list per condition.
#' @return A Scatter plot colored by the expression of the gene list .
#' @examples
#' score_gene_list(CDSC,gene_list,title,searching_mode="match",scat_per_gene=TRUE,normalize="Monocle",n_cores=4)
df_genes_logfc <- function(object,target,signif_genes,colorby,top_n=15,normalize="Monocle",scale="zscore",alter=FALSE,qval_thres=0.05,fc_thres=1,each_cl=TRUE,n_cores=4,batch_formula="~num_genes_expressed"){
    GE_matrix <- data.matrix(object, rownames.force = NA)
    Normilize_Vector <- pData(object)[, 'Size_Factor']
    count_matrix <- as.matrix(sweep(GE_matrix,2,Normilize_Vector,"/"))
    # ==================== Calculate LogFoldChange ===============================
    fc_dataset <- count_matrix[match(signif_genes$gene_short_name,rownames(count_matrix)),]
    rownames(fc_dataset) <- signif_genes$gene_short_name
    # Calculate Mean of each Cluster
    zz<-t(rowsum(t(fc_dataset),group=object[[target]]))
    Mean_exp<- sweep(zz,2,as.vector(table(object[[target]])),"/")
    print("ok")
    # Assign genes to Cluster
    identity <- colnames(Mean_exp)[max.col(Mean_exp)]
    signif_genes["cluster"] <- identity

    # Calculate LogFoldChange
    iter<-0
    LogFoldChange <- apply(fc_dataset,1,function(x){
        iter <<- iter+1
        log2(
        (mean(x[object[[target]]  ==  signif_genes$cluster[iter]]) +1)   /
        (mean(x[object[[target]]  !=  signif_genes$cluster[iter]]) +1)
    )})
    signif_genes["avg_logFC"] <- LogFoldChange
    # all_diff_test_res_condition["identity"] <- levels(object[[target]])[identity]
    all_sig_genes <- signif_genes

    # Subset the genes with qval < threshold
    sig_genes <- subset(signif_genes, qval < qval_thres & avg_logFC >= fc_thres)
    indexes <- sort(sig_genes$qval,index.return=TRUE)$ix
    gnames_pvals <- sig_genes[indexes,c("gene_short_name","pval","qval","avg_logFC","cluster")]
    # --------------------------

    if (each_cl){

        if(get_os()=="windows"){
            n_cores=1
        }
        # if(target=="condition"){
        #   heat_genes <- unlist(mclapply(levels(object[[target]]),function(x){
        #     top_genes <- as.vector(gnames_pvals[levels(object[[target]])[gnames_pvals$Cluster]==x,]$gene_short_name[1:top_n])
        #   },mc.cores=n_cores))
        # }else{
        heat_genes <- unlist(mclapply(levels(object[[target]]),function(x){
            top_genes <- as.vector(gnames_pvals[gnames_pvals$cluster==x,]$gene_short_name[1:top_n])
          },mc.cores=n_cores))
        # }
    }else{
        heat_genes <- as.vector(gnames_pvals$gene_short_name[1:top_n])
    }

    heat_genes<-heat_genes[!is.na(heat_genes)]
    return(list("sig_genes"=heat_genes,"df_pval_genes"=all_sig_genes%>% dplyr::arrange(qval)))

}



