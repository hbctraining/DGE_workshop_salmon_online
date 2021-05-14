# Differential Gene Expression Analysis (bulk RNA-seq Part II)

## Learning Objectives

- Explain and interpret QC on count data using Principal Component Analysis (PCA) and hierarchical clustering
- Implement DESeq2 to obtain a list of significantly different genes
- Perform functional analysis on gene lists with R-based tools

## Installations

[Follow the instructions linked here](../README.md#installation-requirements) to download R and RStudio + Install Packages from CRAN and Bioconductor

## Lessons

### Part 1 (Getting Started)
1. [Workflow (raw data to counts)](../lessons/01a_RNAseq_processing_workflow.md)
1. [Experimental design considerations](../lessons/experimental_planning_considerations.md)
1. [Intro to DGE / setting up DGE analysis](../lessons/01b_DGE_setup_and_overview.md)
     
***

### Part II (QC and setting up for DESeq2)
1. [RNA-seq counts distribution](../lessons/01c_RNAseq_count_distribution.md)
1. [Count normalization](../lessons/02_DGE_count_normalization.md)
1. [Sample-level QC](../lessons/03_DGE_QC_analysis.md) (PCA and hierarchical clustering)
1. [Design formulas](../lessons/04a_design_formulas.md)
1. [Hypothesis testing and multiple test correction](../lessons/05a_hypothesis_testing.md)

***

### Part III (DESeq2)
1. [Description of steps for DESeq2](../lessons/04b_DGE_DESeq2_analysis.md)
1. [Wald test results](../lessons/05b_wald_test_results.md)
1. [Summarizing results and extracting significant gene lists](../lessons/05c_summarizing_results.md)
1. [Visualization](../lessons/06_DGE_visualizing_results.md)
1. [Likelihood Ratio Test results](../lessons/08a_DGE_LRT_results.md)
1. [Time course analysis](../lessons/08b_time_course_analyses.md)

***

### Part IV (Functional Analysis)
1. [Gene annotation](../lessons/genomic_annotation.md)
1. [Functional analysis - over-representation analysis](../lessons/10_FA_over-representation_analysis.md)
1. [Functional analysis - functional class scoring / GSEA](../lessons/11_FA_functional_class_scoring.md)

***
   
[Workflow Summary](../lessons/07_DGE_summarizing_workflow.md)

***

## Building on this workshop
* [Single-cell RNA-seq workshop](https://hbctraining.github.io/scRNA-seq/)
* [RMarkdown](https://hbctraining.github.io/Training-modules/Rmarkdown/)
* [Ggplot2 for functional analysis](https://hbctraining.github.io/Training-modules/Tidyverse_ggplot2/lessons/ggplot2.html)

## Resources
* [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory-behind-deseq2)
* GitHub book on [RNA-seq gene level analysis](http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html)
* [Bioconductor support site](https://support.bioconductor.org/t/deseq2/) (posts tagged with `deseq2`) 
* [Functional analysis visualization](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)

****

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
