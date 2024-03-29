# Mouse microglia regional heterogeneity in single cell resolution

[![DOI](https://zenodo.org/badge/536983198.svg)](https://zenodo.org/badge/latestdoi/536983198)


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
https://doi.org/10.3389/fimmu.2021.639613

Oihane Uriarte Huarte, Dimitrios Kyriakis, Tony Heurtaux, Yolanda Pires-Afonso, Kamil Grzyb, Rashi Halder, Manuel Buttini, Alexander Skupin, Michel Mittelbronn, Alessandro Michelucci

## Overview



### Quality control
(**Folder: QC**)
We filtered out the low quality cells and genes separately in each data set. We defined cells as low-quality, based on three criteria for each cell.  The number of the genes that expressed is more than 200 and 2 median-absolute- deviations (MADs) above the median, the total number of counts is 3 MADs above or below the median and the percentage of counts to mitochondrial genes is 1.5 median-absolute- deviations (MADs) above the median. Cells failing at least two criterion were considered as low quality cells and filtered out from further analysis (Sup. Fig1). Similar to the cell filtering, we filtered out the low quality genes that been expressed in less than 10 cells in the data.


<img src="Plots/QC.png" alt="some text">

We used Seurat v3.1 for the analysis of the gene expression matrix. We merged the gene expression matrices from striatum and midbrain. After the filtering, the count data of 1337 cells and 13446 genes was used for the downstream analysis. To identify the different cell types, we clustered the cells and visualized clusters using a t-distributed stochastic neighbor embedding (t-SNE) plot (Fig. 1 b). The resolution for the Louvain clustering, selected based on the silhouette score performance of different resolutions. This revealed eight main clusters. 

### Cell Identity
(**Folder: Cell Types, DF_Clusters, Garnett**)
In order to recover the cell type identity of each cluster we performed one unsupervised and one supervised workflow.  In the former, we performed a differential expression analysis between the clusters. The genes that popped up as differentially expressed, were searched in literature and we linked them with specific cell types (Fig1 d). In the latter one, we used known cell specific markers (Fig1 e).  These two workflows let us assign each cluster to cell types (Fig1 c,f).  In the folder DF_Clusters there are also barplots for the 3top markers for each cluster.


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
(**Folder: Microglia_Analysis**)
Since we focus here on the microglia, we subset the data and keep only the cells that were identified as microglia. We re calculate the most variable genes only in microglia population and then re-project and cluster these cells. (**Folder: Microglia_Analysis**)
This revealed four different clusters (Fig2, a). The minor cluster (cluster 4) seems to express also oligodendrocyte markers. This also supported by the assignment of the cell identities using garnett.

**&ast;<span style="color:red">The differential expression analysis was based on counts and not on the intergraded data, as suggesterd in https://satijalab.org/seurat/faq.html FAQ4</span>.**


<figure>
    <img src="Plots/Microglia1.png" alt="some text" width=45%>
    <img src="Plots/DF_Micro_1.png" alt="some text" width=45%>
    <figcaption><b>Fig.?:</b> Heatmap of the differentially expressed genes. 
    <b>a)</b> tSNE projection of Microglia 
    <b>b)</b> Differentially expression analysis between the microglia subpopulations.</figcaption>
</figure>

We filtered out the fourth cluster from the further analysis. Then we reprojected with two different ways. The first using RNA counts and one using the Intergrated data (suggested).
Then we performed differential expression analysis between the clusters and the brain regions (Fig2, b). 

### Reproject using Counts (**folder**microglia-heterogeneity/Microglia_Analysis/Reanalysis_counts)
#### Differential Expression Between Clusters
<figure>
    <img src="Plots/DF_Micro_Seurat.png" alt="some text" width=45%>
    <img src="Plots/DF_Micro_Monocle.png" alt="some text" width=45%>
    <figcaption><b>Fig.?:</b> Heatmap of the differentially expressed genes. 
    <b>a)</b> Differentially expression analysis between the microglia subpopulations using Seurat. 
    <b>b)</b> Differentially expression analysis between the microglia subpopulations using Seurat.</figcaption>
</figure>
<br><br><br>
Also, I selected these genes to see the expression of them. Clearly there is a pattern. In monocle some of them poped up as df using the qvalue as criterion of significance, and thery are not the top ones.
<br>
<img src="Plots/Barplot_Homeostatic.png" alt="some text" width=45%>

USING MONOCLE<br>
  status           family         pval         qval gene_short_name cluster avg_logFC
"18"	"OK"	"negbinomial.size"	1.0242427547937e-07	7.64995533080365e-05	"FCRLS"	"2"	0.515234295669162
"33"	"OK"	"negbinomial.size"	8.34449147926016e-06	0.00339949525597496	"LGMN"	"2"	0.280826139509967
"35"	"OK"	"negbinomial.size"	1.13790515893961e-05	0.00437085627336688	"P2RY12"	"2"	0.434463416975197
"538"	"OK"	"negbinomial.size"	0.015847063116479	0.395999844866066	"MARCKSL1"	"3"	0.141552356468414

<br><br>
USING SEURAT<br>
"p_val" "avg_logFC" "pct.1" "pct.2" "p_val_adj" "cluster"   "gene"  "qvalue" <br>

0.232572343774478	0.0631369904109933	0.643	0.838	1	"1"	"P2RY12"
0.48840635534124	0.201733242386338	0.625	0.747	1	"1"	"LGMN"
0.00448532237475512	0.60979247535684	0.717	0.516	1	"2"	"FCRLS"
0.0096385813499733	0.402861609301031	0.906	0.745	1	"2"	"P2RY12"
0.0533954481259887	0.291918466530621	0.811	0.682	1	"2"	"LGMN"
0.00268195851660255	0.892718309072014	0.137	0.019	1	"3"	"MARCKSL1"
0.000286727351104985	0.0885953224050939	0.86	0.475	1	"4"	"FCRLS"
0.00210272546684706	0.0546859483767581	0.94	0.738	1	"4"	"P2RY12"
0.549413314129729	0.106905837504546	0.92	0.65	1	"4"	"LGMN"


#### Differential Expression Between Conditions

<figure>
    <img src="Plots/DF_Cond_Micro_Seurat.png" alt="some text" width=45%>
    <figcaption>Fig.2: Heatmap of the differentially expressed genes between brain regions.</figcaption>
</figure>

#### APOE
<img src="Plots/APOE.png" alt="some text" width=45%>





### Reproject using Intergrated (**folder**microglia-heterogeneity/Microglia_Analysis/Reanalysis)
#### Differential Expression Between Clusters
<figure>
    <img src="Plots/DF_Micro_Seurat2.png" alt="some text" width=45%>
    <img src="Plots/DF_Micro_Monocle2.png" alt="some text" width=45%>
    <figcaption><b>Fig.?:</b> Heatmap of the differentially expressed genes. 
    <b>a)</b> Differentially expression analysis between the microglia subpopulations using Seurat. 
    <b>b)</b> Differentially expression analysis between the microglia subpopulations using Seurat.</figcaption>
</figure>

USING MONOCLE<br>
"95"  "OK"  "negbinomial.size"  0.000446390461852046  0.0631930483073759  "LGMN"  "2" 0.35446083096855
"112" "OK"  "negbinomial.size"  0.000528687200958001  0.0672706468676459  "P2RY12"  "2" 0.446555222869192
"188" "OK"  "negbinomial.size"  0.00239210984441033 0.167196347203629 "FCRLS" "2" 0.334943881552685
"1137"  "OK"  "negbinomial.size"  0.0344099184213774  0.40632487394412  "MARCKSL1"  "3" 0.0662047100564074

<br><br>
USING SEURAT<br>
"p_val" "avg_logFC" "pct.1" "pct.2" "p_val_adj" "cluster"   "gene"  <br>
0.438801091359106 0.220772678403435 0.059 0.04  1 "1" "MARCKSL1"
0.013009652836217 0.373079024728738 0.712 0.715 1 "2" "LGMN"
0.262720171117617 0.410522824906289 0.55  0.577 1 "2" "FCRLS"
0.630789328652911 0.419980206982374 0.762 0.8 1 "2" "P2RY12"
0.0688711980529796  0.439556026513453 0.111 0.03  1 "3" "MARCKSL1"



#### Differential Expression Between Conditions

<figure>
    <img src="Plots/DF_Cond_Micro_Seurat2.png" alt="some text" width=45%>
    <figcaption>Fig.2: Heatmap of the differentially expressed genes between brain regions.</figcaption>
</figure>


#### APOE
<img src="Plots/APOE2.png" alt="some text" width=45%>
