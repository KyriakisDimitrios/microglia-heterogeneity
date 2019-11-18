# Mouse microglia regional heterogeneity in single cell resolution



### RNA-seq analysis

#### R session information

For requirements, please see the file RNA-seq_sessionInfo.txt for specific R packages and their versions.

#### Preprocessing

R code to preprocess single-cell RNA-seq data from scratch can be found in the src directory. Steps (e.g. quality control,normalization, clustering etc.) should be run in the order indicated by the leading counter in the filename (e.g. Outliers_Detectio,.R, Clustering.R )

#### Figures

We additionally provide scripts to reproduce most of the main figures. The R code can be found in the figure_src directory and is split into individual files. Figures can be recreated without the need of preprocessing the data.

### Citation

Please refer to the following research article when using data from this repository:

O. Uriarte, D. Kyriakis, T. Heurtaux, K. Grzyb, R. Halder, M. Buttini, A. Skupin, M. Mittelbronn& A. Michelucci


## Overview





We filtered out the low quality cells and genes separately in each data set. We defined cells as low-quality, based on three criteria for each cell.  The number of the genes that expressed is more than 200 and 2 median-absolute- deviations (MADs) above the median, the total number of counts is 3 MADs above or below the median and the percentage of counts to mitochondrial genes is 1.5 median-absolute- deviations (MADs) above the median. Cells failing at least one criterion were considered as low quality cells and filtered out from further analysis (Sup. Fig1). Similar to the cell filtering, we filtered out the low quality genes that been expressed in less than 10 cells in the data.


<img src="Plots/QC.png" alt="some text">

We used Seurat v3.1 for the analysis of the gene expression matrix. We merged the gene expression matrices from striatum and midbrain. After the filtering, the count data of 1036 cells and 16648 genes was used for the downstream analysis. To identify the different cell types, we clustered the cells and visualized clusters using a t-distributed stochastic neighbor embedding (t-SNE) plot (Fig. 1 b). The resolution for the Louvain clustering, selected based on the silhouette score performance of different resolutions. This revealed eight main clusters. 
In order to recover the cell type identity of each cluster we performed one unsupervised and one supervised workflow.  In the former, we performed a differential expression analysis between the clusters. The genes that popped up as differentially expressed, were searched in literature and we linked them with specific cell types (Fig1 d). In the latter one, we used known cell specific markers (Fig1 e).  These two workflows let us assign each cluster to cell types (Fig1 c,f). 
Since we focus here on the microglia, we subset the data and keep only the cells that were identified as microglia. We re-project and cluster these cells. This revealed four different clusters (Fig2, a). The minor cluster (cluster 4) seems to express also oligodendrocyte markers so we filtered out the four cluster. Then we performed differential expression analysis between the brain regions (Fig2, b).

<figure>
  <img src="Plots/Data.png" alt="some text">
  <img src="Plots/DF_Clusters.png" alt="some text" width=45%>
  <img src="Plots/Heat_Cell_Markes.png" alt="Heat_Cell_Markes" width=45%>
  <img src="Plots/Barplot.png" alt="some text">
  <figcaption>Fig.1: Two-dimensional t-SNE plot of 1036 individual cells at two different regions of mouse brain. Each dot represents a single cell.  
  a) The colors correspond to different regions of the mouse brain. 
  b) t-SNE plot where the colors correspond to eight different clusters. 
  c) The colors correspond to different cell types, which were assigned based on specific cell markers. 
  d) Heat map of the top 15 differentially expressed genes for each cluster. 
  e) Heat map using specific cell type markers.
  f) Bar plot of cell type markers of each cell, colored based on the different cell types.
  g) Bar plot of with the number of cells from the different brain regions per clusters.</figcaption>
</figure>

<figure>
    <img src="Plots/DF_Micro_1.png" alt="some text">
    <img src="Plots/DF_Micro_2.png" alt="some text">
    <img src="Plots/DF_Micro_2_Cond.png" alt="some text">
    <figcaption>Fig.2: Heatmap of the differentially expressed genes. A) Differentially expression analysis between the microglia subpopulations. b) Differentially expression analysis between the striatum and midbrain.</figcaption>
</figure>

