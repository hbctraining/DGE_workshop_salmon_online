---
title: "Set up and overview for gene-level differential expression analysis"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
date: "May 12, 2017"
---

[GEO]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51443 "Gene Expression Omnibus"
[SRA]: https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=SRP031507 "Sequence Read Archive"

Approximate time: 60 minutes

## Learning Objectives 

* Describe the RNA-seq and the differential gene expression analysis workflow
* Explain the experiment and its objectives
* Create a project in R 
* Setup for the analysis of RNA-seq data


# Differential gene expression analysis

Over the past decade, RNA sequencing (RNA-seq) has become an indispensable tool for transcriptome-wide analysis of differential gene expression and differential splicing of mRNAs[1](https://www.nature.com/articles/s41576-019-0150-2). The correct identification of which genes/transcripts are changing in expression between specific conditions is key in our understanding of the biological processes that are affected. 

In this workshop, we will walk you through an **end-to-end gene-level RNA-seq differential expression workflow** using various R packages. We will start with reading in data obtained from Salmon, convert pseudocounts to counts, perform exploratory data analysis for quality assessment and to explore the relationship between samples, perform differential expression analysis, and visually explore the results prior to performing downstream functional analysis.

<p align="center">
<img src="../img/de_workflow_salmon.png" width="400">
</p>


## Review of the dataset

For this workshop, we will be using a publicly available RNA-Seq dataset that is part of a larger study described in [Kenny PJ et al, Cell Rep 2014](http://www.ncbi.nlm.nih.gov/pubmed/25464849). 

The RNA-Seq was performed on HEK293F cells that were either transfected with a MOV10 transgene, or siRNA to knock down Mov10 expression, or non-specific (irrelevant) siRNA. This resulted in 3 conditions **Mov10 oe** (over expression), **Mov10 kd** (knock down) and **Irrelevant kd**, respectively. The number of replicates is as shown below. 

Using these data, we will evaluate transcriptional patterns associated with perturbation of MOV10 expression. Please note that the irrelevant siRNA will be treated as our control condition.

<p align="center">
<img src="../img/dataset.png" width="400">
</p>

***What is the purpose of these datasets? What does Mov10 do?***

The authors are investigating interactions between various genes involved in Fragile X syndrome, a disease in which there is aberrant production of the FMRP protein. 

> **FMRP** is “most commonly found in the brain, is essential for normal cognitive development and female reproductive function. Mutations of this gene can lead to fragile X syndrome, mental retardation, premature ovarian failure, autism, Parkinson's disease, developmental delays and other cognitive deficits.” - from [wikipedia](https://en.wikipedia.org/wiki/FMR1)

> **MOV10**, is a putative RNA helicase that is also associated with **FMRP** in the context of the microRNA pathway. 

**The hypothesis [the paper](http://www.ncbi.nlm.nih.gov/pubmed/25464849) is testing is that FMRP and MOV10 associate and regulate the translation of a subset of RNAs.**

<p align="center">
<img src="../img/mov10-model.png" width="400">
</p>

**Our questions:**
* What patterns of expression can we identify with the loss or gain of MOV10? 
* Are there any genes shared between the two conditions?

### RNA-seq workflow

For this dataset, raw sequence reads were obtained from the [Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/?term=SRP029367). These reads were then processed using the RNA-seq workflow as detailed in the [pre-reading for this workshop](01a_RNAseq_processing_workflow.md). All steps were performed on the command line (Linux/Unix), including a thorough quality control assessment. If you are interested, we have the **MultiQC html report for this dataset [linked here](https://github.com/hbctraining/DGE_workshop_salmon_online/raw/master/data/multiqc_report_rnaseq.html.zip)** for you to peruse.

The directories of output from the mapping/quantification step of the workflow (Salmon) is the data that we will be using. These transcript abundance estimates, often referred to as **'pseudocounts', will be the starting point for our differential gene expression analysis**.



## Setting up

Let's get started by opening up RStudio and setting up a new project for this analysis. 

1. Go to the `File` menu and select `New Project`.
2. In the `New Project` window, choose `New Directory`. Then, choose `New Project`. Name your new directory `DEanalysis` and then "Create the project as subdirectory of:" the Desktop (or location of your choice).
3. The new project should automatically open in RStudio. 

To check whether or not you are in the correct working directory, use `getwd()`. The path `Desktop/DEanalysis` should be returned to you in the console. Within your working directory use the `New folder` button in the bottom right panel to create two new directories:  `meta` and `results`. Remember the key to a good analysis is keeping organized from the start! (*NOTE: we will be downloading our `data` folder`*)

Now we need to grab the files that we will be working with for the analysis. There are two things we need to download. 

1. First we need the **Salmon results for the full dataset**. *Right click on the links below, and choose the "Save link as ..." option to download directly into your project directory*:

* [Salmon data](https://www.dropbox.com/s/oz9yralwbtphw8u/data.zip?dl=1) for the Mov10 full dataset

Once you have the zip file downloaded you will want to decompress it. This will create a `data` directory with sub-directories that correspond to each of the samples in our dataset.

2. Next, we need the **annotation file** which maps our transcript identifiers to gene identifiers. We have created this file for you using the R Bioconductor package [AnnotationHub](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html). For now, we will use it as is but later in the workshop we will spend some time showing you how to create one for yourself. *Right click on the links below, and choose the "Save link as ..." option to download directly into your project directory.*

* [Annotation file](https://github.com/hbctraining/DGE_workshop_salmon/raw/master/data/tx2gene_grch38_ens94.txt)

Finally, go to the `File` menu and select `New File`, then select `R Script`. This should open up a script editor in the top left hand corner. This is where we will be typing and saving all commands required for this analysis. In the script editor type in header lines:

```
## Gene-level differential expression analysis using DESeq2
```

Now save the file as `de_script.R`. When finished your working directory should now look similar to this:


![setup](../img/settingup.png)


### Loading libraries

For this analysis we will be using several R packages, some which have been installed from CRAN and others from Bioconductor. To use these packages (and the functions contained within them), we need to **load the libraries.** Add the following to your script and don't forget to comment liberally!

```r
## Setup
### Bioconductor and CRAN libraries used
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)
```

### Loading data

The main output of Salmon is a `quant.sf` file, and we have one of these for each individual sample in our dataset. An screenshot of the file is displayed below:

<img src="../img/quant_screenshot.png" width="400">

For each transctipy that was assayed in the reference, we have:

1. The transcript identifier
2. The transcript length (in bp)
3. The effective length (described in detail below) 
4. TPM (transcripts per million), which is computed using the effective length
5. The estimated read count ('pseudocount')


> #### What exactly is the effective length?
> The sequence composition of a transcript affects how many reads are sampled from it. While two transcripts might be of identical actual length, depending on the sequence composition we are more likely to generate fragments from one versus the other. The transcript that has a higer likelihood of being sampled, will end up with the larger effective length. The effective lengths is transcript length which has been "corrected" to include factors due to sequence-specific and GC biases.

We will be using the R Bioconductor package `tximport` to prepare the `quant.sf` files for DESeq2. The first thing we need to do is create a variable that contains the paths to each of our `quant.sf` files. Then we will **add names to our quant files which will allow us to easily discriminate between samples in the final output matrix**.


```R
## List all directories containing data  
samples <- list.files(path = "./data", full.names = T, pattern="salmon$")

## Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## Since all quant files have the same name it is useful to have names for each element
names(files) <- str_replace(samples, "./data/", "") %>% 
                str_replace(".salmon", "")
```

Our Salmon index was generated with transcript sequences listed by Ensembl IDs, but `tximport` needs to know **which genes these transcripts came from**. We will use the annotation table that we downloaded to extract transcript to gene information. 

```r
# Load the annotation table for GrCh38
tx2gene <- read.delim("tx2gene_grch38_ens94.txt")

# Take a look at it 
tx2gene %>% View()
```

**`tx2gene`** is a three-column **dataframe linking transcript ID (column 1) to gene ID (column 2)** to gene symbol (column 3). We will take the first two columns as input to `tximport`. The **column names are not relevant, but the column order is (i.e transcript ID must be first).**

Now we are ready to **run `tximport`**. 

```R
?tximport   # let's take a look at the arguments for the tximport function
```

The `tximport()` function imports transcript-level estimates from various external software and summarizes to the gene-level (default) or outputs transcript-level matrices. There are optional arguments to use the abundance estimates as they appear in the `quant.sf` files or to comput alternative values.

For our analysis we **need non-normalized count estimates at the gene-level for performing DESeq2 analysis**. 

Since the gene-level count matrix is a default (`txOut=FALSE`) there is only one additional argument for us to modify **to specify how to obtain our count values**. The options fpr `countsFromAbundance` are as follows:

* `no` (default): This will take the values in TPM and NumReads columns and collapse it down to the gene-level.
* `scaledTPM`: This is taking the TPM scaled up to library size. 
* `lengthScaledTPM`: This is used to generate the count table from the TPM (rather than summarizing the NumReads column). Count values are generated by using the TPM value x featureLength x library size. These represent quantities that are on the same scale as original counts, except no longer correlated with feature length across samples.

> Recall, how the TPM is typically computed:
> 1. Divide the read counts by the length of each gene in kilobases. This gives you reads per kilobase (RPK).
> 2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
> 3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.

```R
# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("tx_id", "ensgene")], countsFromAbundance="lengthScaledTPM")
```

> **An additional argument for `tximport`**: When performing your own analysis you may find that the reference transcriptome file you obtain from Ensembl will have version numbers included on your identifiers (i.e ENSG00000265439.2). This will cause a discrepancy with the tx2gene file since the annotation databases don't usually contain version numbers (i.e ENSG00000265439). To get around this issue you can use the argument `ignoreTxVersion	= TRUE`. The logical value indicates whether to split the tx id on the '.' character to remove version information, for easier matching.

### Viewing data

The `txi` object is a simple list containing matrices of the abundance, counts, length. Another list element 'countsFromAbundance' carries through the character argument used in the tximport call. The length matrix contains the average transcript length for each gene which can be used as an offset for gene-level analysis.

```r
attributes(txi)

$names
[1] "abundance"           "counts"              "length"             
[4] "countsFromAbundance"
```

We will be using the `txi` object as is, for input into DESeq2 but will save it until the next lesson. **For now let's take a look at the count matrix.** You will notice that there are decimal values, so let's round to the nearest whole number and convert it into a dataframe. We wil save it to a variable called `data` that we can play with. 

```r
# Look at the counts
txi$counts %>% View()

# Write the counts to an object
data <- txi$counts %>% 
  round() %>% 
  data.frame()

```

> ### What if I don't have Salmon pseudocounts as input?
> Until recently, the standard approach for RNA-seq analysis had been to map our reads using a splice-aware aligner (i.e STAR) and then use the resulting BAM files as input to counting tools like [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) and [htseq-count](https://htseq.readthedocs.io/en/release_0.11.1/count.html) to obtain our final expression matrix. The field has now moved towards using lightweight alignment tools like Salmon as standard practice. If you are still working with data generated using the older standard approach we have some [materials linked here](https://hbctraining.github.io/DGE_workshop/schedule/1.5-day.html) on **using DESeq2 with a raw count matrix as your starting point**.


### Creating metadata

Of great importance is keeping track of the information about our data. At minimum, we need to at least **have a file which maps our samples to the corresponding sample groups that we are investigating**. We will use the columns headers from the counts matrix as the row names of our metadata file and have single column to identify each sample as "MOV10_overexpression", "MOV10_knockdown", or "control". 

```r
## Create a sampletable/metadata
sampletype <- factor(c(rep("control",3), rep("MOV10_knockdown", 2), rep("MOV10_overexpression", 3)))
meta <- data.frame(sampletype, row.names = colnames(txi$counts))
```

**Now we are all set to start our analysis!**

***

[Next Lesson](01c_RNAseq_count_distribution.md)

***

***
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
