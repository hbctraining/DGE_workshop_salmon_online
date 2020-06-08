---
title: "Exploring DESeq2 results: Wald test"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
date: "June 1, 2020"
---

Approximate time: 60 minutes

## Learning Objectives 

* LFC shrinkage 
* Gene-level filtering?
* Building results tables for comparison of different sample classes
* Summarizing significant differentially expressed genes for each comparison

# Exploring Results (Wald test)

By default DESeq2 uses the Wald test to identify genes that are differentially expressed between two sample classes. Given the factor(s) used in the design formula, and how many factor levels are present, we can extract results for a number of different comparisons. Here, we will walk through how to obtain results from the `dds` object and provide some explanations on how to interpret them.

## Specifying contrasts

In our dataset, we have three sample classes so we can make three possible pairwise comparisons:

1. Control vs. Mov10 overexpression
2. Control vs. Mov10 knockdown
3. Mov10 knockdown vs. Mov10 overexpression

**We are really only interested in #1 and #2 from above**. When we intially created our `dds` object we had provided `~ sampletype` as our design formula, indicating that `sampletype` is our main factor of interest.

To indicate which two sample classes we are interested in comparing, we need to specify **contrasts**. The contrasts are used as input to the DESeq2 `results()` function to extract the desired results. 

Contrasts can be specified in two different ways (with the first method more commonly used):

1. Contrasts can be supplied as a **character vector with exactly three elements**: the name of the factor (of interest) in the design formula, the name of the two factors levels to compare. The factor level given last is the base level for the comparison. The syntax is given below:
	
```r
	# DO NOT RUN!
	contrast <- c("condition", "level_to_compare", "base_level")
	results(dds, contrast = contrast)
```

2. Contrasts can be given as a **list of 2 character vectors**: the names of the fold changes for the level of ineterest, and the names of the fold changes for the base level. These names should match identically to the elements of `resultsNames(object)`. *This method can be useful for combining interaction terms and main effects.*

```r
	# DO NOT RUN!
	resultsNames(dds) # to see what names to use
	contrast <- list(resultsNames(dds)[1], resultsNames(dds)[2])
	results(dds, contrast = contrast)
```

Alternatively, if you **only had two factor levels you could do nothing** and not worry about specifying contrasts (i.e. `results(dds)`). In this case, DESeq2 will choose what your base factor level based on alphabetical order of the levels.

To start, we want to evaluate **expression changes between the MOV10 overexpression samples and the control samples**. As such we will use the first method for specifcying contrasts and create a character vector:

```r
 
## Define contrasts for MOV10 overexpression
contrast_oe <- c("sampletype", "MOV10_overexpression", "control")
```

> ### Does it matter what I choose to be my base level?
> 
> Yes, it does matter. **Deciding what level is the base level will determine how to interpret the fold change that is reported.**  So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in factor level of interest relative to the base level. Thus, if leaving it up to DESeq2 to decide on the contrasts be sure to check that the alphabetical order coincides with the fold change direction you are anticipating.


## The results table

Now that we have our contrast created, we can use it as input to the `results()` function. Let's take a quick look at the help manual for the function:

```r
?results
```
You will see we have the option to provide a wide array of arguments and tweak things from the defaults as needed. As we go through the lesson we will keep coming back to the help doumentation to discuss some arguments that are good to know about.
 
```r
## Extract results for MOV10 overexpression vs control
res_tableOE <- results(dds, contrast=contrast_oe, alpha = 0.05)
```

> **NOTE:** For our analysis, in addition to the `contrast` argument we will also provide a value of 0.05 for the `alpha` argument. We will describe this in more detail when we talk about [gene-level filtering]().

The results table that is returned to us is **a `DESeqResults` object**, which is a simple subclass of DataFrame. In many ways it can be treated like a dataframe (i.e when accessing/subsetting data), however it is important to recognize that there are differences for downstream steps like visualization.

```r
# Check what type of object is returned
class(res_tableOE)
```

Now let's take a look at **what information is stored** in the results:

```r
res_tableOE %>% 
data.frame() %>% 
View()
```

```
log2 fold change (MAP): sampletype MOV10_overexpression vs control 
Wald test p-value: sampletype MOV10_overexpression vs control 
DataFrame with 57914 rows and 6 columns
               		baseMean	log2FoldChange	lfcSE		stat		pvalue		padj
              		<numeric>	<numeric>	<numeric>	<numeric>	<numeric>	<numeric>
ENSG00000000003		3.53E+03	-0.427190489	0.0755347	-5.65604739	1.55E-08	4.47E-07
ENSG00000000005		2.62E+01	0.016159765	0.23735203	0.06584098	9.48E-01	9.74E-01
ENSG00000000419		1.48E+03	0.362663551	0.10761742	3.36995355	7.52E-04	4.91E-03
ENSG00000000457		5.19E+02	0.219135591	0.09768842	2.24476439	2.48E-02	8.21E-02
ENSG00000000460		1.16E+03	-0.261603812	0.07912962	-3.30661411	9.44E-04	5.92E-03
...			...		...		...		...		...		...
```

We have six columns of information reported for each gene (row). We can use the `mcols()` function to extract information on what the values stored in each column represent:

```r
mcols(res_tableOE, use.names=T)
```

* `baseMean`: mean of normalized counts for all samples
* `log2FoldChange`: log2 fold change
* `lfcSE`: standard error
* `stat`: Wald statistic
* `pvalue`: Wald test p-value
* `padj`: BH adjusted p-values
 

## P-values

The p-value is a probability value used to determine whether there is evidence to reject the null hypothesis. **A smaller p-value means that there is stronger evidence in favor of the alternative hypothesis**. However, because we are performing a test for each inidividual gene we need to correct these p-values for multiple testing.

**The `padj` column** in the results table represents the adjusted p-value. The default method for **multiple test correction** in DESeq2 is an implementation of the Benjamini Hochberg false discovery rate (FDR). There are other corrections methods available and can be changed by adding the `pAdjustMethod` argument to the `results()` function.

**The `padj` column is the most important column of the results**. In order to identify a set of genes which are significantly differentially expressed you will want to set a threshold. Typically, `padj` < 0.05 is a good starting point.

### Gene-level filtering

Let's take a closer look at our results table. As we scroll through it, you will notice that for **selected genes there are NA values in the `pvalue` and `padj` columns**. What does this mean?


<p align="center">
<img src="../img/gene_filtering.png" width="700">
</p>

The missing values represent genes that have undergone filtering as part of the `DESeq()` function. Prior to differential expression analysis it is **beneficial to omit genes that have little or no chance of being detected as differentially expressed.** This will increase the power to detect differentially expressed genes. DESeq2 does not physically remove any genes from the original counts matrix, and so all genes will be present in your results table. The genes omitted by DESeq2 meet one of the **three filtering criteria outlined below**:

**1. Genes with zero counts in all samples**

If within a row, all samples have zero counts there is no expression information and therefore these genes are not tested. 

```r
res_tableOE[which(res_tableOE$baseMean == 0),] %>% 
data.frame() %>% 
View()
```

> **The baseMean column for these genes will be zero, and the log2 fold change estimates, p-value and adjusted p-value will all be set to NA.**


**2. Genes with an extreme count outlier**

The `DESeq()` function calculates, for every gene and for every sample, a diagnostic test for outliers called Cook’s distance. Cook’s distance is a measure of how much a single sample is influencing the fitted coefficients for a gene, and a large value of Cook’s distance is intended to indicate an outlier count. Genes which contain a Cook’s distance above a threshold are flagged, however at least 3 replicates are required for flagging, as it is difficult to judge which sample might be an outlier with only 2 replicates. We can turn off this filtering by using the `cooksCutoff` argument in the `results()` function.

```r
res_tableOE[which(res_tableOE$pvalue == NA & res_tableOE$padj == NA),] %>% 
data.frame() %>% 
View()
```

> **If a gene contains a sample with an extreme count outlier then the p-value and adjusted p-value will be set to NA.** 


**3. Genes with a low mean normalized counts**

DESeq2 defines a low mean threshold, that is empirically determined from your data, in which the fraction of significant genes can be increased by reducing the number of genes that are considered for muliple testing. This is based on the notion that genes with very low counts are not likely to see significant differences typically due to high dispersion.

<p align="center">
<img src="../img/indep_filt_scatterplus.png" width="450">
</p>

*Image courtesy of [slideshare presentation](https://www.slideshare.net/joachimjacob/5rna-seqpart5detecting-differentialexpression) from Joachim Jacob, 2014.*

At a user-specified value (`alpha = 0.05`), DESeq2 evaluates the change in the number of significant genes as it filters out increasingly bigger portions of genes based on their mean counts, as shown in the figure above. The point at which the number of significant genes reaches a peak is the low mean threshold that is used to filter genes that undergo multiple testing. There is also an argument to turn off the filtering off by setting `independentFiltering = F`.

```r
res_tableOE[which(res_tableOE$pvalue != NA & res_tableOE$padj == NA),] %>% 
data.frame() %>% 
View()
```

> **If a gene is filtered by independent filtering, then only the adjusted p-value will be set to NA. 

> **NOTE:** DESeq2 will perform the filtering outlined above by default; however other DE tools, such as EdgeR will not.  Filtering is a necessary step, even if you are using limma-voom and/or edgeR's quasi-likelihood methods. Be sure to follow pre-filtering steps when using other tools, as outlined in their user guides found on Bioconductor as they generally perform much better. 

## Fold change

**The order of the names in the contrast determines the direction of fold change that is reported.** The name provided in the second element is the level that is used as baseline. So for example, if we observe a log2 fold change of -2 this would mean the gene expression is lower in Mov10_oe relative to the control. However, these estimates do not account for the large dispersion we observe with low read counts. To avoid this, the **log2 fold changes calculated by the model need to be adjusted**. 

>
> **NOTE:** The Wald test can also be used with **continuous variables**. If the variable of interest provided in the design formula is continuous-valued, then the reported log2 fold change is per unit of change of that variable.
>


### More accurate LFC estimates

To generate more accurate log2 foldchange estimates, DESeq2 allows for the **shrinkage of the LFC estimates toward zero** when the information for a gene is low, which could include:

- Low counts
- High dispersion values

As with the shrinkage of dispersion estimates, LFC shrinkage uses **information from all genes** to generate more accurate estimates. Specifically, the distribution of LFC estimates for all genes is used (as a prior) to shrink the LFC estimates of genes with little information or high dispersion toward more likely (lower) LFC estimates. 

<img src="../img/deseq2_shrunken_lfc.png" width="500">

*Illustration taken from the [DESeq2 paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8).*

For example, in the figure above, the green gene and purple gene have the same mean values for the two sample groups (C57BL/6J and DBA/2J), but the green gene has little variation while the purple gene has high levels of variation. For the green gene with low variation, the **unshrunken LFC estimate** (vertex of the green **solid line**) is very similar to the shrunken LFC estimate (vertex of the green dotted line), but the LFC estimates for the purple gene are quite different due to the high dispersion. So even though two genes can have similar normalized count values, they can have differing degrees of LFC shrinkage. Notice the **LFC estimates are shrunken toward the prior (black solid line)**.

In the most recent versions of DESeq2, the shrinkage of LFC estimates is **not performed by default**. This means that the log2 foldchanges would be the same as those calculated by:

```r
log2 (normalized_counts_group1 / normalized_counts_group2)
```

To generate the shrunken log2 fold change estimates, you have to run an additional step on your results object (that we will create below) with the function `lfcShrink()`.

```r
## Save the unshrunken results to compare
res_tableOE_unshrunken <- res_tableOE

# Apply fold change shrinkage
res_tableOE <- lfcShrink(dds, contrast=contrast_oe, res=res_tableOE)
```

> **NOTE: Shrinking the log2 fold changes will not change the total number of genes that are identified as significantly differentially expressed.** The shrinkage of fold change is to help with downstream assessment of results. For example, if you wanted to subset your significant genes based on fold change for further evaluation, you may want to use shruken values. Additionally, for functional analysis tools such as GSEA which require fold change values as input you would want to provide shrunken values.

### Adding a fold change threshold: 
With large significant gene lists it can be hard to extract meaningful biological relevance. To help increase stringency, one can also **add a fold change threshold**. Use shrunken values here!

For e.g., we can create a new threshold `lfc.cutoff` and set it to 0.58 (remember that we are working with log2 fold changes so this translates to an actual fold change of 1.5).

`lfc.cutoff <- 0.58`

`sigOE <- res_tableOE_tb %>% filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)`

> ### An alternative approach to add the fold change threshold:
> The `results()` function has an option to add a fold change threshold using the `lfcThrehsold` argument. This method is more statistically motivated, and is recommended when you want a more confident set of genes based on a certain fold-change. It actually performs a statistical test against the desired threshold, by performing a two-tailed test for log2 fold changes greater than the absolute value specified. The user can change the alternative hypothesis using `altHypothesis` and perform two one-tailed tests as well. **This is a more conservative approach, so expect to retrieve a much smaller set of genes!**
>
> Test this out using our data:
> 
> `results(dds, contrast = contrast_oe, alpha = 0.05, lfcThreshold = 0.58)`
>
> **How do the results differ? How many significant genes do we get using this approach?**


## Visualizing results with an MA plot

A plot that can be useful to exploring our results is the MA plot. The MA plot shows the mean of the normalized counts versus the log2 foldchanges for all genes tested. The genes that are significantly DE are colored to be easily identified. This is also a great way to illustrate the effect of LFC shrinkage. The DESeq2 package offers a simple function to generate an MA plot. 

**Let's start with the unshrunken results:**

```r
plotMA(res_tableOE_unshrunken, ylim=c(-2,2))
```

<img src="../img/maplot_unshrunken.png" width="600">

**And now the shrunken results:**

```r
plotMA(res_tableOE, ylim=c(-2,2))
```

<img src="../img/MA_plot.png" width="600">

In addition to the comparison described above, this plot allows us to evaluate the magnitude of fold changes and how they are distributed relative to mean expression. Generally, we would expect to see significant genes across the full range of expression levels. 


***

**Excerise**

**MOV10 Differential Expression Analysis: Control versus Knockdown**

Now that we have results for the overexpression results, let's do the same for the **Control vs. Knockdown samples**. Use contrasts in the `results()` to extract a results table and store that to a variable called `res_tableKD`.  

```r
## Define contrasts, extract results table and shrink log2 fold changes
contrast_kd <-  c("sampletype", "MOV10_knockdown", "control")

res_tableKD <- results(dds, contrast=contrast_kd, alpha = 0.05)

res_tableKD <- lfcShrink(dds, contrast=contrast_kd, res=res_tableKD)
```

Take a quick peek at the results table containing Wald test statistics for the Control-Knockdown comparison we are interested in and make sure that format is similar to what we observed with the OE.

***

## Summarizing results

To summarize the results table, a handy function in DESeq2 is `summary()`. Confusingly it has the same name as the function used to inspect data frames. This function when called with a DESeq results table as input, will summarize the results using the alpha threshold: FDR < 0.05 (padj/FDR is used even though the output says `p-value < 0.05`). Let's start with the OE vs control results:

```r
## Summarize results
summary(res_tableOE, alpha = 0.05)
```

In addition to the number of genes up- and down-regulated at the default threshold, **the function also reports the number of genes that were tested (genes with non-zero total read count), and the number of genes not included in multiple test correction due to a low mean count**.


## Extracting significant differentially expressed genes

Let's first create variables that contain our threshold criteria:

```r
### Set thresholds
padj.cutoff <- 0.05
```

We can easily subset the results table to only include those that are significant using the `filter()` function, but first we will convert the results table into a tibble:

```r
res_tableOE_tb <- res_tableOE %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
```

Now we can subset that table to only keep the significant genes using our pre-defined thresholds:

```r
sigOE <- res_tableOE_tb %>%
        filter(padj < padj.cutoff)
```

```r
sigOE
```

***

**Exercise**

Do the same with KD
Using the same p-adjusted threshold as above (`padj.cutoff < 0.05`), subset `res_tableKD` to report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.

```r

res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()
  
sigKD <- res_tableKD_tb %>%
        filter(padj < padj.cutoff)
```

**How many genes are differentially expressed in the Knockdown compared to Control?** 
```r
sigKD
``` 
***


Now that we have extracted the significant results, we are ready for visualization!



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

*Some materials and hands-on activities were adapted from [RNA-seq workflow](http://www.bioconductor.org/help/workflows/rnaseqGene/#de) on the Bioconductor website*

***
