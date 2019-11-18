#' Outlier Detection.
#' @author Dimitrios Kyriakis
#' @export
#' @param df: Gene expression matrix
#' @param outlier_detector: Method to detect the outliers "MAD"/"MEAN" (Default="MAD").
#'
#' @return Gene expression Matrix without Outliers
#' @examples
#' Outlier_detection(df,outlier_detector="MAD")
Outlier_detection <- function(object,outlier_detector="MAD",min.features=200,min.cells=5) {
    cat(green("Outlier Detection Started"))
    nCounts <- object$nCount_RNA
    nFeatures <- object$nFeature_RNA
    percent.mit <- object$percent.mito
    percent.rb <- object$percent.rb

    if (tolower(outlier_detector)=="mad"){
        print("MAD outlier Detection")
        false_list <- rep(TRUE,length(nFeatures))
        res_out_nFeat <- outliers_mad(nFeatures,threshold = 2)
        out_nFeat<-false_list
        out_nFeat[res_out_nFeat$outliers_pos] <-FALSE


        res_out_nCounts <- outliers_mad(nCounts,threshold = 3)
        out_nCounts <- res_out_nCounts$outliers_pos
        out_nCount <- false_list
        out_nCount[out_nCounts]<-FALSE

        res_out_percent.mit <- outliers_mad(percent.mit,threshold = 1.5)
        out_percent.mit <- percent.mit[res_out_percent.mit$outliers_pos]> mean(percent.mit)
        out_percent.mits <- res_out_percent.mit$outliers_pos[out_percent.mit]
        out_percent.mit <- false_list
        out_percent.mit[out_percent.mits] <-FALSE


        res_out_nCounts <- res_out_nCounts$limit
        res_out_nFeat <- c(min.features,res_out_nFeat$limit[2])
        res_out_percent.mit <-res_out_percent.mit$limit[2]

        object <- object[,out_nFeat & out_nCount & out_percent.mit ]
        oultliers_index <- out_nFeat & out_nCount & out_percent.mit

    }else if (tolower(outlier_detector)=="mean"){
        # We filter cells that have unique feature counts over +2*sd or less than 200
        # We filter cells that have >2*sd% mitochondrial counts
        print("MEAN outlier Detection")
        outliers_nFeature_RNA <- (nFeatures > min.features) & (nFeatures < round(mean(nFeatures)+2*sd(nFeatures)))
        outliers_nCount_RNA <- (nCounts < round(mean(nCounts)+2*sd(nCounts)))
        outliers_percent.mt <- (percent.mit < mean(percent.mit)+2*sd(percent.mit))

        object <- subset(x = object, subset = outliers_nFeature_RNA & outliers_nCount_RNA & outliers_percent.mt )
        oultliers_index <- outliers_nFeature_RNA & outliers_nCount_RNA & outliers_percent.mt

        res_out_nCounts <- c(min.features,round(mean(nFeatures)+2*sd(nFeatures)))
        res_out_nFeat <- c(round(mean(nCounts)-2*sd(nCounts)),round(mean(nCounts)+2*sd(nCounts)))
        res_out_percent.mit <- mean(percent.mit)+2*sd(percent.mit)

    }else{
        df_filtered <- object
        oultliers_index <- NULL
        return(list("object_filtered"=df_filtered,"oultliers_index"=oultliers_index))
    }
    # ======================= OUTLIER PLOT HISTOFGRAMMS ==============================
    pdf(paste(Sys.Date(),"_Outliers_hist_",object$stim,".pdf",sep=""),width=7,height=3)
    par(mfrow=c(1,3))
    hist(nCounts, col="lightblue",main="Histogram for nCounts"); abline(v=res_out_nCounts,col="red", lwd=3, lty=2)
    hist(nFeatures, col="lightblue",main="Histogram for nFeatures"); abline(v=res_out_nFeat,col="red", lwd=3, lty=2)
    hist(percent.mit, col="lightblue",main="Histogram for percent.mit");  abline(v=res_out_percent.mit,col="red", lwd=3, lty=2)
    dev.off()
    # ----------------------------------------------------------------------------------
    print("Dimensions after Removing Outliers")
    print(dim(object))
    cat(green("Outlier Detection Finished\n"))
    # print("Outlier Detection Finished")
    return(list("object_filtered"=object,"oultliers_index"=oultliers_index))
}

