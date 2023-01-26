# Differential gene expression workshop using Salmon counts

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R/) | 4-session online workshop (~8 hours of trainer-led time)|

### Description

This repository has teaching materials for a hands-on **Introduction to differential gene expression (DGE) analysis** workshop. The workshop will lead participants through performing a differential gene expression analysis workflow on RNA-seq count data using R/RStudio. Working knowledge of R is required or completion of the [Introduction to R workshop](https://hbctraining.github.io/Intro-to-R/). 

**Note for Trainers:** Please note that the schedule linked below assumes that learners will spend between 3-4 hours on reading through, and completing exercises from selected lessons between classes. The online component of the workshop focuses on more exercises and discussion/Q & A.

### Learning Objectives

- QC on count data using Principal Component Analysis (PCA) and hierarchical clustering
- Using DESeq2 to obtain a list of significantly different genes
- Visualizing expression patterns of differentially expressed genes
- Performing functional analysis on gene lists with R-based tools

### Lessons
* [Workshop schedule (trainer-led learning)](schedule/)
* [Self-learning](schedule/links-to-lessons.md)

### Installation Requirements

Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) (version 4.0.0 or above)
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)
 
Note:  When installing the following packages, if you are asked to select (a/s/n) or (y/n), please select “a” or "y" as applicable.

(1) Install the below packages on your laptop from CRAN. You DO NOT have to go to the CRAN webpage; you can use the following function to install them one by one:


```r
install.packages("insert_package_name_in_quotations")
install.packages("insert_package_name_in_quotations")
& so on ...
```

Note that these package names are case sensitive!

```r
BiocManager
tidyverse
RColorBrewer
pheatmap
ggrepel
cowplot
```

(2) Install the below packages from Bioconductor. Load BiocManager, then run BiocManager's `install()` function 12 times for the 12 packages:

```r
library(BiocManager)
install("insert_first_package_name_in_quotations")
install("insert_second_package_name_in_quotations")
& so on ...
```

Note that these package names are case sensitive!

```r
DESeq2
clusterProfiler
DOSE
org.Hs.eg.db
pathview
DEGreport
tximport
AnnotationHub
ensembldb
```

> **NOTE:** The library used for the annotations associated with genes (here we are using `org.Hs.eg.db`) will change based on organism (e.g. if studying mouse, would need to install and load `org.Mm.eg.db`). The list of different organism packages are given [here](https://github.com/hbctraining/Training-modules/raw/master/DGE-functional-analysis/img/available_annotations.png).

(3) Finally, please check that all the packages were installed successfully by loading them **one at a time** using the code below:  

```r
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(cowplot)
library(clusterProfiler)
library(DEGreport)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(tximport)
library(AnnotationHub)
library(ensembldb)
```

(4) Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```

---

### Citation

To cite material from this course in your publications, please use:

> Meeta Mistry, Mary Piper, Jihe Liu, & Radhika Khetani. (2021, May 24). hbctraining/DGE_workshop_salmon_online: Differential Gene Expression Workshop Lessons from HCBC (first release). Zenodo. https://doi.org/10.5281/zenodo.4783481

A lot of time and effort went into the preparation of these materials. Citations help us understand the needs of the community, gain recognition for our work, and attract further funding to support our teaching activities. Thank you for citing this material if it helped you in your data analysis.

---

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
