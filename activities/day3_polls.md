1. Please choose the option that is **NOT TRUE**.

	1. Generally, the gene dispersion estimates should decrease with increasing mean expression values.
	1. Generally, the gene dispersion estimates should cluster around the line of best fit.
	1. Gene dispersion plots are a good way to examine whether the data are a good fit to the DESeq2 model.
	1. **A good gene dispersion plot indicates a higher likelihood of detecting a lot of differentially expressed genes.**
	1. A worrisome gene dispersion plot could be a result of outlier samples or contamination.
	1. A worrisome gene dispersion plot indicates that we should be more skeptical of our significant DE genes.

2. In order to set up your null hypothesis, you need to first observe all of the data points in an experiment.

	1. True
	1. **False**

3. The generalized linear model (GLM) is fit for each individual gene, which means we are conducting thousands of independent tests for a given experiment. This inflates our false positive rate for DE genes.

	1. **True**
	1. False

4. Which of the following statements about gene-level filtering is **NOT TRUE** ?

	1. If a gene has zero expression in all of the samples it is not tested for differential expression.
	1. For independent filtering that is applied in DESeq2, the low mean threshold is empirically determined from your data.
	1. **Gene-level filtering will reduce the number of genes being tested and therefore decrease the total number of differentially expressed genes that are identified.**
	1. A Cook's distance is computed for each gene in each sample to help identify genes with an extreme outlier count.

5. Which of the following statements about LFC shrinkage is true?

	1. LFC shrinkage uses information from the significant genes to generate more accurate estimates.
	1. Shrinking the log fold changes will reduce the total number of genes identified as significant at padj < 0.05.
	1. **Shrinkage of the LFC estimates is useful when the information for a gene is low, which includes low mean expression or high dispersion**.
	1. LFC shrinkage is applied by default in DESeq2. 


6. To identify significant differentially expressed genes, you will need to set an adjusted p-value threshold and a fold change threshold. 

	1. True
	1. **False**
	
7. When using a heatmap to observe the differences in patterns of gene expression between samplegroups it is more informative to plot the normalized expression values, rather than the scaled Z-scores for each gene.

	1. True
	1. **False**

