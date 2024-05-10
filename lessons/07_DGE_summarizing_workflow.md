---
title: "Summary of DGE workflow"
author: "Mary Piper"
date: "June 8, 2017"
---

Approximate time: 15 minutes

## Learning Objectives 

* Identify the R commands needed to run a complete differential expression analysis using DESeq2

## Summary of differential expression analysis workflow

We have detailed the various steps in a differential expression analysis workflow, providing theory with example code. To provide a more succinct reference for the code needed to run a DGE analysis, we have summarized the steps in an analysis below:

1. Obtaining gene-level counts from Salmon using tximport

	```r
	# Run tximport
	txi <- tximport(files, 
			type="salmon", 
			tx2gene=t2g, 
			countsFromAbundance = "lengthScaledTPM")
	
	# "files" is a vector wherein each element is the path to the salmon quant.sf file, and each element is named with the name of the sample.
	# "t2g" is a 2 column data frame which contains transcript IDs mapped to geneIDs (in that order)
	```

2. Creating the dds object:
		
	```r
	# Check that the row names of the metadata equal the column names of the **raw counts** data
	all(colnames(txi$counts) == rownames(metadata))
	
	# Create DESeq2Dataset object
	dds <- DESeqDataSetFromTximport(txi, 
					colData = metadata, 
					design = ~ condition)
	```
	
3. Exploratory data analysis (PCA & hierarchical clustering) - identifying outliers and sources of variation in the data:
	
	```r
	# Transform counts for data visualization
	rld <- rlog(dds, 
		    blind=TRUE)
	
	# Plot PCA 
	plotPCA(rld, 
		intgroup="condition")
	
	# Extract the rlog matrix from the object and compute pairwise correlation values
	rld_mat <- assay(rld)
	rld_cor <- cor(rld_mat)
	
	# Plot heatmap
	pheatmap(rld_cor, 
		 annotation = metadata)
	```
	
4. Run DESeq2:

	```r
		# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis in include other sources of variation using "dds <- DESeqDataSetFromTximport(txi, colData = metadata, design = ~ covariate + condition)"

	# Run DESeq2 differential expression analysis
	dds <- DESeq(dds)

		# **Optional step** - Output normalized counts to save as a file to access outside RStudio using "normalized_counts <- counts(dds, normalized=TRUE)"
	```
	
5. Check the fit of the dispersion estimates:
	
	```r
	# Plot dispersion estimates
	plotDispEsts(dds)
	``` 

6. Create contrasts to perform Wald testing on the shrunken log2 foldchanges between specific conditions:

	```r
	# Specify contrast for comparison of interest
	contrast <- c("condition", "level_to_compare", "base_level")
	
	# Output results of Wald test for contrast
	res <- results(dds, 
		       contrast = contrast, 
		       alpha = 0.05)
	
	# Shrink the log2 fold changes to be more accurate
	res <- lfcShrink(dds, 
			 coef = "sampletype_group1_vs_group2", 
			 type = "apeglm")	 
         # The coef will be dependent on what your contrast was. and should be identical to what is stored in resultsNames()
	```

7. Output significant results:

	```r
	# Set thresholds
	padj.cutoff < - 0.05
	
	# Turn the results object into a tibble for use with tidyverse functions
	res_tbl <- res %>%
                  data.frame() %>%
                  rownames_to_column(var="gene") %>% 
                  as_tibble()
	
	# Subset the significant results
	sig_res <- dplyr::filter(res_tbl, 
			  padj < padj.cutoff)
	```

8. Visualize results: volcano plots, heatmaps, normalized counts plots of top genes, etc.

9. Perform analysis to extract functional significance of results: GO or KEGG enrichment, GSEA, etc.

10. Make sure to output the versions of all tools used in the DE analysis:

	```r
	sessionInfo()
	```
	
	For better reproducibility, it can help to create **RMarkdown reports**, which save all code, results, and visualizations as nicely formatted html reports. We have a very basic example of a report [linked here](https://www.dropbox.com/s/4bq0chxze6dogba/workshop-example.html?dl=0). To create these reports we have [additional materials](https://hbctraining.github.io/Training-modules/Rmarkdown/) available.
