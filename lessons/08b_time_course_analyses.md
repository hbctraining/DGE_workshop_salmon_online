---
title: "Time course analysis with DESeq2"
author: "Meeta Mistry and Mary Piper"
date: "June 17, 2020"
---

Approximate time: 20 minutes

## Learning Objectives 

* Discuss time course analyses with DESeq2

## Time course analyses

The LRT test can be especially helpful when performing time course analyses. We can use the LRT to explore whether there are any significant differences in treatment effect between any of the timepoints. 

For have an experiment looking at the effect of treatment over time on mice of two different genotypes. We could use a design formula for our 'full model' that would include the major sources of variation in our data: `genotype`, `treatment`, `time`, and our main condition of interest, which is the difference in the effect of treatment over time (`treatment:time`).

```r
full_model <- ~ genotype + treatment + time + treatment:time
```

To perform the LRT test, we can determine all genes that have significant differences in expression between treatments across time points by giving the 'reduced model' without the `treatment:time` term:

```r
reduced_model <- ~ genotype + treatment + time
```

Then, we could run our test by using the following code:

```r
dds <- DESeqDataSetFromMatrix(countData = raw_counts, colData = metadata, design = ~ genotype + treatment + time + treatment:time)

dds_lrt_time <- DESeq(dds, test="LRT", reduced = ~ genotype + treatment + time)
```
This analysis will not return genes where the treatment effect does not change over time, even though the genes may be differentially expressed between groups at a particular time point, as shown in the figure below:

<img src="../img/lrt_time_nodiff.png" width="200">

The significant DE genes will represent those genes that have differences in the effect of treatment over time, an example is displayed in the figure below:

<img src="../img/lrt_time_yesdiff.png" width="200">


Once we have our results, we can determine the significant genes using a threshold of `padj` < 0.05 and return the normalized counts for those genes. Then we could perform clustering to identify genes that change over time in a way meaningful to us:

```r
clusters <- degPatterns(cluster_rlog, metadata = meta, time="time", col="treatment")
```

You can extract the groups of genes associated with the patterns of interest similar to the actions performed previously, then move on to functional analysis for each of the gene groups of interest.

## Resources

We have covered the inner workings of DESeq2 in a fair amount of detail such that when using this package you have a good understanding of what is going on under the hood. For more information on topics covered, we encourage you to take a look at the following resources:

* [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory-behind-deseq2)
* GitHub book on [RNA-seq gene level analysis](http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html)
* [Bioconductor support site](https://support.bioconductor.org/t/deseq2/) (posts tagged with `deseq2`) 

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*
