# Microglia-heterogeneity

Mouse microglia regional heterogeneity in single cell resolution

## RNA-seq analysis
R session information
For requirements, please see the file RNA-seq_sessionInfo.txt for specific R packages and their versions.

## Preprocessing
R code to preprocess single-cell RNA-seq data from scratch can be found in the src directory. Steps (e.g. quality control, filtering, normalization, etc.) should be run in the order indicated by the leading counter in the filename (e.g. 00_filtering.R before 01_normalisation.R)

### Figures
We additionally provide scripts to reproduce most of the main figures. The R code can be found in the figure_src directory and is split into individual files. Figures can be recreated without the need of preprocessing the data.


### Figures
We provide scripts to reproduce figures from the ATAC-seq data analysis of the manuscript. The R code can be found in the figure_src directory.

## Citation
O. Uriarte, D. Kyriakis, T. Heurtaux, K. Grzyb, R. Halder, M. Buttini, A. Skupin, M. Mittelbronn& A. Michelucci