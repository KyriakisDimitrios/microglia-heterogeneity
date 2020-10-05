#' Load Files and Create Matrix.
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param list_of_files: a list of files to be load.
#' @param iter_qc: Iteration.
#' @param data_10x: Type of data (Default 10x=FALSE).
#'
#' @return Gene expression Matrix
#' @examples load_files(list_of_files,iter_qc,data_10x=FALSE)
load_files <- function(file,data_10x){
    if(data_10x){
        print(file)
        barcode.path <- paste0(file, "barcodes.tsv")
        features.path <- paste0(file, "features.tsv")
        matrix.path <- paste0(file, "matrix.mtx")
        mat <- as.matrix(readMM(file = matrix.path))
        feature.names = read.delim(features.path, header = FALSE, stringsAsFactors = FALSE)
        barcode.names = read.delim(barcode.path, header = FALSE, stringsAsFactors = FALSE)
        colnames(mat) = barcode.names$V1
        rownames(mat) = feature.names$V2
        object <- mat[,order(colSums(as.matrix(mat)),decreasing=T)]
    }else{
        object <- read.csv(file,header=T,row.names=1, sep='\t')
        object <- object[,order(colSums(object),decreasing=T)]
    }
    # = Change the cases to UPPER and - to dots
    rows_nms <- str_replace_all( toupper(rownames(object)),"-",".")
    rownames(object) <- rows_nms
    return(object)
}


#' remove_genes
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: Gene expression matrix.
#'
#' @return Gene expression Matrix
#' @examples remove_genes(object)
remove_genes <- function(object,remove_mt=TRUE,remove_rb=TRUE){
    # ========= Mitochondria Percentage ==========
    Mit.index <- grep(pattern = "^MT-|^MT\\.", x = rownames(object), value = FALSE)
    percent.Mit <- Matrix::colSums(object[Mit.index, ])/Matrix::colSums(object)
    # --------------------------------------------

    # ========= Ribosomal Percentage ==========
    RIB.index <- grep(pattern = "^RPL|^RPS", x = rownames(object), value = FALSE)
    percent.RIB <- Matrix::colSums(object[RIB.index, ])/Matrix::colSums(object)
    # -----------------------------------------

    if(remove_mt){
        object <- object[-Mit.index, ]
    }

    if(remove_rb){
        object <- object[-RIB.index, ]
    }

    return(object)
}


#' Metrics Calculation Mit,RB,ERCC
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param object: Gene expression matrix.
#'
#' @return Gene expression Matrix
#' @examples metrics_calc(object)
metrics_calc <- function(object,remove_mt=TRUE,remove_rb=TRUE){
    # ============= ERCC Percentage ==============
    ERCC.index <- grep(pattern = "^ERCC", x = rownames(object), value = FALSE)
    percent.ERCC <- Matrix::colSums(object[ERCC.index, ])
    object <- object[-ERCC.index, ]
    # --------------------------------------------

    # ========= Mitochondria Percentage ==========
    Mit.index <- grep(pattern = "^MT-|^MT\\.", x = rownames(object), value = FALSE)
    percent.Mit <- Matrix::colSums(object[Mit.index, ])/Matrix::colSums(object)
    # --------------------------------------------


    if(remove_mt){
        object <- object[-Mit.index, ]
    }

    # ========= Ribosomal Percentage ==========
    RIB.index <- grep(pattern = "^RPL|^RPS", x = rownames(object), value = FALSE)
    percent.RIB <- Matrix::colSums(object[RIB.index, ])/Matrix::colSums(object)
    # -----------------------------------------

    if(remove_rb){
        object <- object[-RIB.index, ]
    }

    return(list("object"=object,"percent.ERCC"=percent.ERCC,"percent.mito"=percent.Mit,"percent.rb"=percent.RIB))
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







#' Calculates how many Principal Components contain N\% of the variation of the data.


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


#' Read Files
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param list_of_files: elist of files to be load.
#' @param condition_names: call rate.
#' @param remove_mt: Remove mitochondrial genes (Default=TRUE).
#' @param gene_list_remove : Remove list of genes.
#' @param imputation: Impute missing gene expressions with SAVER (Default=FALSE).
#' @param na_threshold: Genes with number of missing values more than threshold deleted.
#' @param n_cores: number of cores to use for Parallel Quality Control.
#' @param elbow: Filtering based on the PCA variance explain elbow.
#' @param data_10x: Type of data (Default 10x=FALSE).
#' @param outlier_detector: Method to detect the outliers "MAD"/"MEAN" (Default="MAD").
#'
#' @return Monocle Cell Data Set
#' @return Gene expression Matrix
#' @examples create_cds(list_of_files,condition_names,data_10x=FALSE,outlier_detector="MAD")
read_file <- function(iter_qc,list_of_files,condition_names,min.features=200,min.cells=5,remove_mt=TRUE,remove_rb=TRUE,outlier_detector="MAD",data_10x=FALSE,elbow=FALSE,imputation=FALSE,tool="Seurat"){
    init_num <- iter_qc
    file <- list_of_files[iter_qc]
    # ============== Load File ===============
    object <- load_files(file = file,data_10x=data_10x)
    # ----------------------------------------

    #============== ELBOW ====================
    if(elbow==TRUE){
        object <- elbow_calc(object,condition_names,iter_qc)
    }
    # ----------------------------------------



    # ================================= Outlier Detection ===================================
    # debugonce(Advanced_Outlier_detection)
    Outlier_object <- Advanced_Outlier_detection(object,condition= condition_names[init_num])
    object <- Outlier_object$object_filtered
    oultliers_index <- Outlier_object$oultliers_index
    # ---------------------------------------------------------------------------------------

    # #============== Metrics Calculation  ====================
    metrics_output <- metrics_calc(object=object,remove_mt=remove_mt,remove_rb=remove_rb)
    object <- metrics_output$object
    percent.mito <- metrics_output$percent.mito
    percent.rb <-  metrics_output$percent.rb
    # -------------------------------------------------------


    Seurat <- CreateSeuratObject(counts = object, project = condition_names[init_num],
                                 min.cells = min.cells,
                                 min.features = min.features,
                                 meta.data = data.frame(percent.mito = percent.mito,percent.rb=percent.rb))
    Seurat$stim <- condition_names[init_num]




    # ================= Standard pre-processing workflow (https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html)
    # QC and selecting cells for further analysis
    # Seurat[["percent.mt"]] <- PercentageFeatureSet(object = Seurat, pattern = "^MT-")
    # pdf(paste(Sys.Date(),"_QC_",Seurat$stim,"before.pdf",sep=""))
    # print(VlnPlot(object = Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.rb"), ncol = 4))
    # dev.off()
    #
    #
    # # # ================================= Outlier Detection ===================================
    # # # debugonce(Outlier_detection)
    # Outlier_object <- Outlier_detection(Seurat,outlier_detector=outlier_detector)
    # Seurat <- Outlier_object$object_filtered
    # oultliers_index <- Outlier_object$oultliers_index
    # # # ---------------------------------------------------------------------------------------

    # =================. Impute missing gene expressions with SAVER =========================
    if (imputation ==TRUE) {
        cat(green("\nImputation with SAVER\n"))
        allcells <- as.matrix(Seurat@assays$RNA@counts)
        library(SAVER)
        allcells[is.na(allcells)] <- 0
        cortex.saver <- saver(allcells, ncores = n_cores)
        allcells<-as.matrix(cortex.saver$estimate)
        allcells[is.na(allcells)] <- 0

        Seurat_imputed <- CreateSeuratObject(counts = allcells, project = condition_names[init_num],
                                     min.cells = min.cells,
                                     min.features = min.features,
                                     meta.data = data.frame(percent.mito = Seurat$percent.mito,percent.rb=Seurat$percent.rb))
        Seurat_imputed$stim <- condition_names[init_num]
        Seurat <- Seurat_imputed

    }
    # --------------------------------------------------------------------------------------

#
#     pdf(paste(Sys.Date(),"_QC_",Seurat$stim,"after.pdf",sep=""))
#     print(VlnPlot(object = Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito","percent.rb"), ncol = 4))
#     dev.off()
    # ----------------------------------------------

    # Seurat <- SCTransform(object = Seurat)
    # ====================== Normalization ====================
    Seurat <- NormalizeData(object = Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
    # ----------------------------------------------


    # ====== Identification of highly variable features (feature selection)
    Seurat <- FindVariableFeatures(object = Seurat, selection.method = "vst", nfeatures = 5000)

    #Seurat <- ScaleData(object = Seurat, vars.to.regress = c("nUMI", "percent.mito"), display.progress = FALSE)

    # Identify the 10 most highly variable genes
    top10 <- head(x = VariableFeatures(object = Seurat), 10)

    # plot variable features with and without labels
    plot1 <- VariableFeaturePlot(object = Seurat)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    # CombinePlots(plots = list(plot1, plot2))
    pdf(paste(Sys.Date(),"_top10_var_feat_",Seurat$stim,".pdf"))
    print(plot2)
    dev.off()
    # ----------------------------------------------


    # Seurat[["percent.OPC"]]     <- PercentageFeatureSet(Seurat, features = c("OLIG1","OLIG2","PDGFRA"))
    # Seurat[["percent.Oligo"]]   <- PercentageFeatureSet(Seurat, features = c("MAG","MOG","PLP1","ERMN","ENPP6"))
    # Seurat[["percent.Astro"]]   <- PercentageFeatureSet(Seurat, features = c("AQP4","ALDOC","MLC1","SLC1A3","GJA1","NDRG2"))
    # Seurat[["percent.Micro"]]   <- PercentageFeatureSet(Seurat, features = c("TMEM119","TREM2","CCR5","LAPTM5","CTSS","AIF1","C1QA"))
    # Seurat[["percent.Neurons"]] <- PercentageFeatureSet(Seurat, features = c("GAD1","GAD2"))
    # Seurat[["percent.NSC"]]     <- PercentageFeatureSet(Seurat, features = c("TOP2A","BIRC5","CDK1"))


    # to merge sets to get all the relevant genes (to be filtered out later again!)
    # BUT HERE ONLY FOR THE 800 MOST EXPRESSING BARCODES!!! TOO EMPTY BEADS CONFUSES MONOCLE!!!
    # tm1 <- Org1[,c(1:number_of_cells)]
    tm1 <- Seurat
    print(dim(tm1))
    # tm_list[[init_num]]  <- tm1
    # init_num <- init_num + 1
    return(tm1)
}



#' Create a newCellDataSet  Monocle object.
#' @author Dimitrios Kyriakis
#' @export
#'
#' @param list_of_files: elist of files to be load.
#' @param condition_names: call rate.
#' @param remove_mt: Remove mitochondrial genes (Default=TRUE).
#' @param gene_list_remove : Remove list of genes.
#' @param imputation: Impute missing gene expressions with SAVER (Default=FALSE).
#' @param na_threshold: Genes with number of missing values more than threshold deleted.
#' @param n_cores: number of cores to use for Parallel Quality Control.
#' @param elbow: Filtering based on the PCA variance explain elbow.
#' @param data_10x: Type of data (Default 10x=FALSE).
#' @param outlier_detector: Method to detect the outliers "MAD"/"MEAN" (Default="MAD").
#'
#' @return Monocle Cell Data Set
#' @return Gene expression Matrix
#' @examples create_cds(list_of_files,condition_names,data_10x=FALSE,outlier_detector="MAD")
create_cds <- function(list_of_files,condition_names,
    min.features=200,min.cells=5,
    remove_mt=TRUE,remove_rb=TRUE,
    outlier_detector="MAD",data_10x=FALSE,elbow=FALSE,imputation=FALSE,
    tool="Seurat",
    n_cores=2){

    # ============================== Checking the Parameters ====================================
    tic()
    cat(cyan("Selected Parameteres :\n")%+% paste("\tRemove Mitochondrial Genes : ",remove_mt,
                                                  "\n\tRemove Ribosomal genes : ",remove_rb,
                                                  "\n\tImputation : ",imputation,
                                                  "\n\t10xData : ",data_10x,
                                                  "\n\tElbow : ",elbow,
                                                  "\n\tOutlier Detector : ",outlier_detector) )
    cat(cyan("\n\nLoading Data : ")  %+% as.character(Sys.time()) %+% "\n" )
    tm_list <- list()
    init_num <- 1

    print("QC for every file")

    if (get_os()=="windows"){
        cat(yellow("Warning :")  %+% "Your os is Windows. The Parallel file reading is implemented only for UNIX/MAC\n")
        n_cores=1
    }
    # debugonce(read_file)
    tm_list <- mclapply(FUN=read_file,c(1:length(list_of_files)),
                                     list_of_files=list_of_files,
                                     condition_names=condition_names,
                                     data_10x=data_10x,
                                     elbow=elbow,remove_mt=remove_mt,remove_rb=remove_rb,
                                    imputation=imputation,
                                     outlier_detector=outlier_detector,
                                      min.features=min.features,
                                      min.cells=min.cells,
                                      mc.cores=n_cores)

    # =========================== Perform integration ==============================
    Seurat.anchors <- FindIntegrationAnchors(object.list = tm_list, dims = 1:20)
    Seurat.combined <- IntegrateData(anchorset = Seurat.anchors, dims = 1:20)
    # doi: https://doi.org/10.1101/576827
    # int.features <- SelectIntegrationFeatures(object.list = tm_list, nfeatures = 3000)
    # tm_list <- PrepSCTIntegration(object.list = tm_list, anchor.features = int.features,
    #                                     verbose = FALSE)
    # int.anchors <- FindIntegrationAnchors(object.list = tm_list, normalization.method = "SCT",
    #                                            anchor.features = int.features, verbose = FALSE)
    # Seurat.combined <- IntegrateData(anchorset = int.anchors, normalization.method = "SCT",
    #                                      verbose = FALSE)


    # =========================== Perform an integrated analysis ==========================
    DefaultAssay(object = Seurat.combined) <- "integrated"
    Seurat.combined$condition <- Idents(object = Seurat.combined)
    # -------------------------------------------------------------------------------------

    # ======================================== CELL Cycle ==========================================
    s.genes <- cc.genes$s.genes
    g2m.genes <- cc.genes$g2m.genes
    if(!g2m.genes%in%rownames(Seurat.combined@assays$integrated@data)){
        if(!s.genes%in%rownames(Seurat.combined@assays$integrated@data)){
        }
    }
    Seurat.combined <- CellCycleScoring(Seurat.combined, s.features = s.genes, g2m.features = g2m.genes)
    # -----------------------------------------------------------------------------------------------------


    object <- Seurat.combined
    if(tolower(tool)=="monocle"){
       object <- seurat_to_monocle(object)
    }
    return(list("Combined"=object,"Data_List"=tm_list))
}
# ----------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------












#' Advanced_Outlier_detection
#' @author Dimitrios Kyriakis
#' @export
#' @param df: Gene expression matrix
#' @param outlier_detector: Method to detect the outliers "MAD"/"MEAN" (Default="MAD").
#'
#' @return Gene expression Matrix without Outliers
#' @examples
#' Advanced_Outlier_detection(df,outlier_detector="MAD")
Advanced_Outlier_detection <- function(object,filtering=TRUE,min.features=200,min.cells=5,condition="") {
    cat(green("Outlier Detection Started"))
    #object <- object1[,!colSums(object1)==0]
    # Count how many genes has non zero value
    nFeatures <- colSums(object != 0)
    #object <- object[,nFeatures>=min.features]
    #object <- object[rowSums(object!=0)>=min.cells,]
    object<- object[rowSums(object!=0)>=min.cells,nFeatures>=min.features]
    nCounts <- colSums(object)
    nFeatures <- colSums(object != 0)
    percent.mit   <-  colSums(object[grep("^MT-|^MT\\.",rownames(object)),])/nCounts
    percent.rb   <-  colSums(object[grep("^RPL|^RPS",rownames(object)),])/nCounts

    df_qc <- data.frame(Cond=rep(condition,dim(object)[2]),
                        nFeatures = nFeatures,
                        Total_MRNA=nCounts,
                        Percent_Mit=percent.mit)
    # nCounts <- object$nCount_RNA
    # nFeatures <- object$nFeature_RNA
    # percent.mit <- object$percent.mito
    # percent.rb <- object$percent.rb

    if(filtering){
        print("MAD outlier Detection")
        false_list <- rep(TRUE,length(nFeatures))
        res_out_nFeat <- outliers_mad(nFeatures,threshold = 2)
        out_nFeat<-false_list
        out_nFeat[res_out_nFeat$outliers_pos] <-FALSE

        # For high count filtering, not exceed the exprected double rate (0.9 per 1000 cells)
        res_out_nCounts <- outliers_mad(nCounts,threshold = 2)
        out_nCounts <- res_out_nCounts$outliers_pos
        out_nCount <- false_list
        out_nCount[out_nCounts]<-FALSE

        if(data_10x){
            doublerate <- round(dim(object)[2]*0.9/1000)
        }

        res_out_percent.mit <- outliers_mad(percent.mit,threshold = 1.5)
        out_percent.mit <- percent.mit[res_out_percent.mit$outliers_pos]> mean(percent.mit)
        out_percent.mits <- res_out_percent.mit$outliers_pos[out_percent.mit]
        out_percent.mit <- false_list
        out_percent.mit[out_percent.mits] <-FALSE

        QC_Matrix <- data.frame("nCounts"=out_nCount,"nFeatures"=out_nFeat,"percent.mit"=out_percent.mit)

        # If a cell drop more than 2 test assinged as outlier (Luecken 2019: Current Best practices in scRNAseq)
        outlier_index <- rep(FALSE,length(nFeatures))
        QC_Matrix$Passed_Tests <- rowSums(QC_Matrix)
        outlier_index[QC_Matrix$Passed_Tests>=2] <- TRUE

        df_filtered <- object[,outlier_index]

        if(res_out_nCounts$limit[1]<0){
            res_out_nCounts$limit[1]<-0
        }
        res_out_nCounts <- res_out_nCounts$limit
        res_out_nFeat <- c(min.features,res_out_nFeat$limit[2])
        res_out_percent.mit <-res_out_percent.mit$limit[2]

    }else{
        df_filtered <- object
        oultliers_index <- NULL
        return(list("object_filtered"=df_filtered,"oultliers_index"=oultliers_index))
    }

    # ======================= OUTLIER PLOT HISTOFGRAMMS ==============================
    pdf(paste(Sys.Date(),"_Outliers_hist_",condition,".pdf",sep=""),width=7,height=3)
    par(mfrow=c(1,3))
    hist(nCounts, col="lightblue",main="Histogram for nCounts"); abline(v=res_out_nCounts,col="red", lwd=3, lty=2)
    hist(nFeatures, col="lightblue",main="Histogram for nFeatures"); abline(v=res_out_nFeat,col="red", lwd=3, lty=2)
    hist(percent.mit, col="lightblue",main="Histogram for percent.mit");  abline(v=res_out_percent.mit,col="red", lwd=3, lty=2)
    dev.off()
    # ----------------------------------------------------------------------------------
    print("Dimensions after Filtering out Low Quality Cells")
    print(dim(df_filtered))

    cat(green("Outlier Detection Finished\n"))

    # ======================= VIOLIN PLOTS ==============================
    # debugonce(Vln_QC)
    Vln_QC(df_qc,condition=condition,title="before",outliers=outlier_index)
    nCounts <- colSums(df_filtered)
    nFeatures <- colSums(df_filtered != 0)
    percent.mit   <-  colSums(df_filtered[grep("^MT-|^MT\\.",rownames(df_filtered)),])/nCounts
    percent.rb   <-  colSums(df_filtered[grep("^RPL|^RPS",rownames(df_filtered)),])/nCounts
    df_qc <- data.frame(Cond=rep(condition,dim(df_filtered)[2]),
                        nFeatures = nFeatures,
                        Total_MRNA=nCounts,
                        Percent_Mit=percent.mit)
    Vln_QC(df_qc,condition=condition,title="after")
    # ------------------------------------------------------------------

    # print("Outlier Detection Finished")
    return(list("object_filtered"=df_filtered,"oultliers_index"=QC_Matrix))
}



#' Elbow Calculation
#' @author Dimitrios Kyriakis
#'
#' @param object: Gene expression matrix.
#' @param condition_names: Condition names.
#' @param iter_qc: Iteration.
#'
#' @return Gene expression Matrix
#' @examples elbow_calc(object,condition_names,iter_qc)
elbow_calc <- function(object,condition_names,iter_qc){
    require(ecp)
    #====================== ELBOW ==========================================
    pdf(paste(Sys.Date(),"QC_Kneeplot_Zeros",condition_names[iter_qc],".pdf",sep="_"))
    # Knee plot
    knee_data <- object
    # plot(colSums(is.na(knee_data)),main="Number of missing genes per cell",ylab="Number of missing genes",xlab="Cells")
    knee_data[is.na(knee_data)] <- 0
    plot(colSums(knee_data == 0),main=condition_names[iter_qc],ylab="Number of non expressed genes",xlab="Cells")
    # print("No memory error")
    # ============== Elbow Analysis ============= #
    
    y<- cumsum(sort(colSums(knee_data),TRUE))
    x<-c(1:length(y))
    y1<- as.matrix(y)
    y2<- as.matrix(y[500:length(y)])
    
    real_index <- e.divisive(diff(y1),k=1,min.size=2)
    real_index <- real_index$considered.last
    print(paste("The real elbow value is ",real_index))
    
    # pdf("knee.pdf")
    plot(x, y, pch=19,xlab="Cells",ylab="Cumulative Total mRNA")
    if (real_index<500){
      y2<- as.matrix(y[real_index:length(y)])
      indicies <- e.divisive(diff(y2),k=1,min.size=2)$considered.last +500
      print(paste("The elbow value adjusted to be between 500 and 1000. The new value is ",indicies))
      points(x[indicies], y[indicies], pch=19, col='red')
    }else{
      indicies <- real_index
    }
    
    points(x[real_index], y[real_index], pch=19, col='lightblue')
    abline(v=500,col="red")
    abline(v=1000,col="red")
    # dev.off()
    if(indicies>1000& real_index<1000){
      object<- object[,1:700]
    }else{
      object<- object[,1:indicies]
    }
    
    # ---------------------------------------------
    dev.off()
    return(object)
}





#' Vln_QC
#' @author Dimitrios Kyriakis
#' @export
Vln_QC<-function(df_qc,condition,title,outliers=NULL) {
  library(gridExtra)
  if(is.null(outliers)){
    outliers <- 1
  }else{
    outliers <- as.numeric(!outliers)+1
  }
  p1 <- ggplot(df_qc,aes(factor(Cond),nFeatures,fill=Cond))+
    geom_violin(width=0.8,show.legend=FALSE)+
    ggtitle("nFeatures") +
    xlab("")+
    ylab("Expression Level")+
    geom_jitter(width = 0.3,size=1,show.legend=FALSE,colour=outliers)
  p2 <- ggplot(df_qc,aes(factor(Cond),Total_MRNA,fill=Cond),width=1)+
    geom_violin(width=0.8,scale="count",show.legend=FALSE,adjust=1/2)+
    ggtitle("nCounts") +
    xlab("")+
    ylab("Expression Level")+
    geom_jitter(width = 0.3,size=1,show.legend=FALSE,colour=outliers)


  p3 <- ggplot(df_qc,aes(factor(Cond),Percent_Mit,fill=Cond),width=1)+
    geom_violin(width=0.8,show.legend=FALSE)+
    ggtitle("Percent Mit") +
    xlab("")+
    ylab("Expression Level")+
    geom_jitter(width = 0.3,size=1,show.legend=FALSE,colour=outliers)

  pdf(paste(Sys.Date(),"QC",condition,title,".pdf",sep="_"))
  grid.arrange(p1, p2, p3,nrow=1)
  dev.off()
}

