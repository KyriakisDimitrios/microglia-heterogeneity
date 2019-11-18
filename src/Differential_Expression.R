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
df_genes <- function(object,target,top_n=15,logfc.threshold=0.25,n_cores=4){
    tool <- object_identifier(object)
    if(tolower(tool)=="monocle"){
        all_diff_test_res_condition <- differentialGeneTest(object, fullModelFormulaStr=paste("~",target,sep=""),cores=n_cores, reducedModelFormulaStr ="~num_genes_expressed")
        # == Reduce Matrix (remove not significant genes)
        signif_genes <-all_diff_test_res_condition
        monocle_helper <- df_genes_logfc(object,target,signif_genes=signif_genes,top_n=15,qval_thres=0.05,fc_thres=logfc.threshold,each_cl=TRUE,n_cores=4)
        pbmc.markers <- monocle_helper$df_pval_genes
        colnames(pbmc.markers)[colnames(pbmc.markers)=="gene_short_name"] <- "gene"
        top10_genes <- monocle_helper$sig_genes

    }else{
        Idents(object) <- object[[target]]
        # find markers for every cluster compared to all remaining cells, report only the positive ones
        pbmc.markers <- FindAllMarkers(object = object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = logfc.threshold)
        top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = top_n, wt = avg_logFC)
        top10_genes<- top10$gene
    }



    annotated_heat(object=object,
                  row_annotation=c(1),
                  gene_list=top10_genes,
                  gene_list_name="DF_genes",
                  title="DF_genes",
                  ordering=target)

    write.table(pbmc.markers, file = paste0(Sys.Date(),"_TO_EXP_each_",target,".tsv"),row.names=FALSE, na="", sep="\t")
    return(list("degs"=pbmc.markers,"top_markers"=top10_genes))
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
