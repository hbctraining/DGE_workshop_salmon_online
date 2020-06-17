---
title: "DGE analysis using LRT in DESeq2"
author: "Meeta Mistry and Mary Piper"
date: "June 14, 2017"
---

Approximate time: 60 minutes

## Learning Objectives 

* Extract results using the LRT and compare to Wald test
* Export results to file


## Exploring results from the Likelihood ratio test (LRT)

DESeq2 also offers the Likelihood Ratio Test as an alternative **when evaluating expression change across more than two levels**. Genes which are identified as significant, are those that are changing in expression in any direction across the different factor levels.

Generally, this test will result in a larger number of genes than the individual pair-wise comparisons. While the LRT is a test of significance for differences of any level(s) of the factor, one should not expect it to be exactly equal to the union of sets of genes using Wald tests (although we do expect a high degree of overlap).

## The `results()` table

To extract the results from our `dds_lrt` object we can us the same `results()` function we had used with the Wald test. _There is no need for contrasts since we are not making a pair-wise comparison._

```r
# Extract results for LRT
res_LRT <- results(dds_lrt)
```

Let's take a look at the results table:

```r
# View results for LRT
res_LRT  
```

```r
log2 fold change (MLE): sampletype MOV10 overexpression vs control 
LRT p-value: '~ sampletype' vs '~ 1' 
DataFrame with 57761 rows and 6 columns
                        baseMean     log2FoldChange              lfcSE             stat               pvalue                 padj
                       <numeric>          <numeric>          <numeric>        <numeric>            <numeric>            <numeric>
ENSG00000000003 3525.88347786355 -0.438245423329571 0.0774607246185232 40.4611749305021 1.63669402960044e-09 3.14070461117016e-08
ENSG00000000005 26.2489043110535 0.0292079869376247  0.441128912409409 1.61898836146221    0.445083140923522     0.58866891597654
ENSG00000000419 1478.25124052691  0.383635036119846  0.113760957175207 11.3410110249776  0.00344612277761083   0.0122924772964227
ENSG00000000457  518.42202383345  0.228970583496456  0.102331174090148 14.6313920603898 0.000665018279181725  0.00304543241149833
ENSG00000000460 1159.77613645835 -0.269138203013482 0.0814992499897986 25.0394477225533 3.65386933066256e-06 3.23415706764646e-05
...                          ...                ...                ...              ...                  ...                  
```

The results table output looks similar to the Wald test results, with identical columns to what we observed previously. 

### Why are fold changes reported for an LRT test?

For analyses using the likelihood ratio test, the p-values are determined solely by the difference in deviance between the full and reduced model formula. **A single log2 fold change is printed in the results table for consistency with other results table outputs, but is not associated with the actual test.**

**Columns relevant to the LRT test:**

* `baseMean`: mean of normalized counts for all samples
* `stat`: the difference in deviance between the reduced model and the full model
* `pvalue`: the stat value is compared to a chi-squared distribution to generate a pvalue
* `padj`: BH adjusted p-values

**Additional columns:**

* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error

> **NOTE:** Printed at the top of the the results table are the two sample groups used to generate the log2 fold change values that we observe in the results table. This can be controlled using the `name` argument; the value provided to name must be an element of resultsNames(dds).

## Identifying significant genes

When filtering significant genes from the LRT we threshold only the `padj` column. _How many genes are significant at `padj < 0.05`?_

```r
# Create a tibble for LRT results
res_LRT_tb <- res_LRT %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset to return genes with padj < 0.05
sigLRT_genes <- res_LRT_tb %>% 
  filter(padj < padj.cutoff)

# Get number of significant genes
nrow(sigLRT_genes)

# Compare to numbers we had from Wald test
nrow(sigOE)
nrow(sigKD)

```

The number of significant genes observed from the LRT is quite high. This list includes genes that can be changing in any direction across the three factor levels (control, KO, overexpression). To reduce the number of significant genes, we can increase the stringency of our FDR threshold (`padj.cutoff`).

***

**Exercise**

1. Compare the resulting gene list from the LRT test to the gene lists from the Wald test comparisons.
    1. How many of the `sigLRT_genes` overlap with the significant genes in `sigOE`?
    1. How many of the `sigLRT_genes` overlap with the significant genes in `sigKD`?

***

### Identifying gene clusters exhibiting particular patterns across samples

Often we are interested in genes that have particular patterns across the sample groups (levels) of our condition. For example, with the MOV10 dataset, we may be interested in genes that exhibit the lowest expression for the `Mov10_KD` and highest expression for the `Mov10_OE` sample groups (i.e. KD < CTL < OE). To identify genes associated with these patterns we can use a clustering tool, `degPatterns` from the 'DEGreport' package, that groups the genes based on their changes in expression across sample groups.

First we will subset our rlog transformed normalized counts to contain only the differentially expressed genes (padj < 0.05).

```r
# Subset results for faster cluster finding (for classroom demo purposes)
clustering_sig_genes <- sigLRT_genes %>%
  arrange(padj) %>%
  head(n=1000)


# Obtain rlog values for those significant genes
cluster_rlog <- rld_mat[clustering_sig_genes$gene, ]
```

Then we can use the `degPatterns` function from the 'DEGreport' package to determine sets of genes that exhibit similar expression patterns across sample groups. The `degPatterns` tool uses a hierarchical clustering approach based on pair-wise correlations, then cuts the hierarchical tree to generate groups of genes with similar expression profiles. The tool cuts the tree in a way to optimize the diversity of the clusters, such that the variability inter-cluster > the variability intra-cluster.

```r
# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = meta, time = "sampletype", col=NULL)
```

Let's explore the output:

```r
# What type of data structure is the `clusters` output?
class(clusters)

# Let's see what is stored in the `df` component
head(clusters$df)
```

While we don't see any clusters with the pattern we are looking for (KD < CTL < OE), we do see a lot of genes that don't change much between control and knockdown sample groups, but increase drastically with the overexpression group (Group 1). 

<img src="../img/degReport_clusters2.png" width="600">

Let's explore the set of genes in Group 1 in more detail:

```r
# Extract the Group 1 genes
cluster_groups <- clusters$df
group1 <- clusters$df %>%
          filter(cluster == 1)
```

After extracting a group of genes, we can perform functional analysis to explore associated functions. We can repeat this extraction and functional analysis for any of the groups of interest.



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *Materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*

