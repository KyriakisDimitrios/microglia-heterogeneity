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
    color_cond  <- c(brewer.pal(5,"Set1"),"black","gray","magenta4","seagreen4")

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
              scale_fill_manual(color_cond)+
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


#' reduce_dim
#' @author Dimitrios Kyriakis
#' @export
reduce_dim<-function(object,project,resolution=NA){
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
        object <- ScaleData(object = object, verbose = FALSE)
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
    optimal_output<- optimal_clusters(object,k.max=10,save=TRUE,resolution=FALSE)
    object <- optimal_output$object
    opt_num <- optimal_output$opt_num
    sil_scor <- optimal_output$sil_scor
    # ---------------------------------------------------------------------------

    plot_cells(object,target="condition",leg_pos="right",save=TRUE,ncol=1)
    plot_cells(object,target="Cluster",leg_pos="right",save=TRUE,ncol=1)

    return(list("Combined"=object,"Num_dim"=num_dim))
}
