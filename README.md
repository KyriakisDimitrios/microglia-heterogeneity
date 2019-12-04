# Mouse microglia regional heterogeneity in single cell resolution



### RNA-seq analysis

#### R session information

For requirements, please see the file RNA-seq_sessionInfo.txt for specific R packages and their versions.
Something to add

#### Preprocessing

R code to preprocess single-cell RNA-seq data from scratch can be found in the src directory. Steps (e.g. quality control,normalization, clustering etc.) should be run in the order indicated by the leading counter in the filename (e.g. Outliers_Detectio,.R, Clustering.R )

#### Figures

We additionally provide scripts to reproduce most of the main figures.
The R code can be found in the figure_src directory and is split into individual files.
Figures can be recreated without the need of preprocessing the data. **NOT READY YET**

### Citation

Please refer to the following research article when using data from this repository:

O. Uriarte


## Overview



### Quality control

We filtered out the low quality cells and genes separately in each data set. We defined cells as low-quality, based on three criteria for each cell.  The number of the genes that expressed is more than 200 and 2 median-absolute- deviations (MADs) above the median, the total number of counts is 3 MADs above or below the median and the percentage of counts to mitochondrial genes is 1.5 median-absolute- deviations (MADs) above the median. Cells failing at least two criterion were considered as low quality cells and filtered out from further analysis (Sup. Fig1). Similar to the cell filtering, we filtered out the low quality genes that been expressed in less than 10 cells in the data.


<img src="Plots/QC.png" alt="some text">

We used Seurat v3.1 for the analysis of the gene expression matrix. We merged the gene expression matrices from striatum and midbrain. After the filtering, the count data of 1337 cells and 13446 genes was used for the downstream analysis. To identify the different cell types, we clustered the cells and visualized clusters using a t-distributed stochastic neighbor embedding (t-SNE) plot (Fig. 1 b). The resolution for the Louvain clustering, selected based on the silhouette score performance of different resolutions. This revealed eight main clusters. 

### Cell Identity

In order to recover the cell type identity of each cluster we performed one unsupervised and one supervised workflow.  In the former, we performed a differential expression analysis between the clusters. The genes that popped up as differentially expressed, were searched in literature and we linked them with specific cell types (Fig1 d). In the latter one, we used known cell specific markers (Fig1 e).  These two workflows let us assign each cluster to cell types (Fig1 c,f). 


<figure>
  <img src="Plots/Data.png" alt="some text">
  <img src="Plots/DF_Clusters.png" alt="some text" width=45%>
  <img src="Plots/Heat_Cell_Markes.png" alt="Heat_Cell_Markes" width=45%>
  <img src="Plots/Barplot.png" alt="some text" width=45%>
  <img src="Plots/Bar_Plots.png" alt="some text" width=45%>
  <figcaption><b>Fig.1:</b> Two-dimensional t-SNE plot of 1036 individual cells at two different regions of mouse brain. Each dot represents a single cell.  
  <b>a)</b> The colors correspond to different regions of the mouse brain. 
  <b>b)</b> t-SNE plot where the colors correspond to eight different clusters. 
  <b>c)</b> The colors correspond to different cell types, which were assigned based on specific cell markers. 
  <b>d)</b> Heat map of the top 15 differentially expressed genes for each cluster. 
  <b>e)</b> Heat map using specific cell type markers.
  <b>f)</b> Bar plot of cell type markers of each cell, colored based on the different cell types.
  <b>g)</b> Bar plot of with the number of cells from the different brain regions per clusters.</figcaption>
</figure>

<figure>
    <img src="Plots/Garnett.png" alt="some text" width=45%>
    <figcaption>t-SNE plot where the colors correspond to garnett cell identity assignment using as cell type markers (<b>Oligodendrocytes:</b> Mog, Mag, Plp1, Mbp, Ermn, Enpp6, Pdgfra, Olig1, Olig2, Cspg4 , <b>Astrocytes:</b> Slc1a2, Aqp4, Gja1, <b>Microglia:</b> Tmem119, P2ry13, Aif1, Cts3, Cx3Xr1, <b>Neurons:</b> Syn1, Gad1, Gad2, Th, Meg3, Slc17a6
, <b>Endothelial:</b> Pecam1, Dcn, Igfbpl1, <b>Ependymal:</b> Igf2, Ccdc153).</figcaption>
</figure>


### Microglia Subpopulations
Since we focus here on the microglia, we subset the data and keep only the cells that were identified as microglia.
We re calculate the most variable genes only in microglia population and then re-project and cluster these cells. 
This revealed four different clusters (Fig2, a). The minor cluster (cluster 4) seems to express also oligodendrocyte markers. This also supported by the assignment of the cell identities using garnett.

**&ast;<span style="color:red">The differential expression analysis was based on counts and not on the intergraded data, as suggesterd in https://satijalab.org/seurat/faq.html FAQ4</span>.**


<figure>
    <img src="Plots/Microglia1.png" alt="some text" width=45%>
    <img src="Plots/DF_Micro_1.png" alt="some text" width=45%>
    <figcaption><b>Fig.?:</b> Heatmap of the differentially expressed genes. 
    <b>a)</b> tSNE projection of Microglia 
    <b>b)</b> Differentially expression analysis between the microglia subpopulations.</figcaption>
</figure>

We filtered out the fourth cluster from the further analysis. Then we performed differential expression analysis between the clusters and the brain regions (Fig2, b). 


#### Differential Expression Between Clusters
<figure>
    <img src="Plots/DF_Micro_Seurat.png" alt="some text" width=45%>
    <img src="Plots/DF_Micro_Monocle.png" alt="some text" width=45%>
    <figcaption><b>Fig.?:</b> Heatmap of the differentially expressed genes. 
    <b>a)</b> Differentially expression analysis between the microglia subpopulations using Seurat order on the average logfc. 
    <b>b)</b> Differentially expression analysis between the microglia subpopulations using Seurat order on the average qvalue.</figcaption>
</figure>
<br><br><br>
Also, I selected these genes to see the expression of them. Clearly there is a pattern. In monocle some of them poped up as df using the qvalue as criterion of significance, and thery are not the top ones.
<br>
<img src="Plots/Barplot_Homeostatic.png" alt="some text" width=45%>

USING MONOCLE<br>
       "status"    "family"    "pval"  "qval"  "gene_short_name"<br>
"LGMN"  "OK"    "negbinomial.size"  8.31768090541528e-06    0.00319517770780881 "LGMN"<br>
"P2RY12"    "OK"    "negbinomial.size"  5.28026686435796e-06    0.0022901028384288  "P2RY12"<br>
"FCRLS" "OK"    "negbinomial.size"  9.85748011699133e-08    6.97546421962886e-05    "FCRLS"<br>
"CSF1"  "OK"    "negbinomial.size"  2.09939470193106e-06    0.00104542080620234 "CSF1"<br>
<br><br>
USING SEURAT<br>
0.00617135162185148 0.404239157791022   0.711   0.709   1   "2" "LGMN"  0.129730086103788<br>



#### Differential Expression Between Conditions

<figure>
    <img src="Plots/DF_Cond_Micro_Seurat.png" alt="some text" width=45%>
    <figcaption>Fig.2: Heatmap of the differentially expressed genes between brain regions order by logfoldChange.</figcaption>
</figure>



