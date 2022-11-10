# 2022HPC
Hypopharyngeal carcinoma (HPC) is a relatively rare malignancy and frequently diagnosed at advanced stage with poor prognosis. Here we provided a comprehensive and unique resource revealing the landscape of HPC TME at single-cell resolution. Through systematic analyses based on scRNA-seq and bulk RNA-seq data from clinical samples, we also uncovered functional gene modules from malignant tumor cells associated with prognosis, and established the relationship between non-malignant subtypes’ composition and patients’ responses to the combined therapy.

## Environment

Linux Server with Ubuntu 20.04.3, R version 3.6.3, Python version 3.7.6

## Install software

### Install R package Seurat v2.3.4

```
remotes::install_version("Seurat", version = "3.2.2")
```

### Install R package DoubletDecon

```
if(!require(devtools)){
  install.packages("devtools") # If not already installed
}
devtools::install_github('EDePasquale/DoubletDecon')
```

### Install R package inferCNV

```
if (!requireNamespace("BiocManager", quietly = TRUE))
     install.packages("BiocManager")
BiocManager::install(version = "devel")
BiocManager::install("infercnv")
```

### Install R package NMF v0.23

```
intall.pacakges('NMF')
```

### Install R package Monocle v2.14

```
source("http://bioconductor.org/biocLite.R") 
biocLite("monocle")	
```

### Install software Scrublet

```
Download from https://github.com/swolock/scrublet
```

### Insatll software CellPhoneDB v2.0.6

```
Download from https://github.com/Teichlab/cellphonedb
```

### Online software Metascape

```
https://metascape.org/gp/index.html	
```

### Online software CIBERSORTx

```
https://cibersortx.stanford.edu/
```





