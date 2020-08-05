1. The following statements are describing the characteristics of RNA-seq count data. Select the statement that is **NOT** correct.     
    1. Count data has large dynamic range, with no upper limit for expression.  
    1. For genes with low mean expression we observe a range of variance values.  
    1. **For genes with high mean expression, the variance is less than the mean.** 
    1. Count data follows a negative binomial distribution.

2. Increasing the number of replicates enables more precise estimates of group means, and increases our statistical power to call differentially expressed genes correctly.  
    1. **True**  
    1. False
   
3. Which of the following factors does NOT need to be taken into consideration when comparing expression of a given gene between samples? 
  
    1. Sequencing depth  
    1. **Gene length**  
    1. RNA composition

4. Which of the following normalization methods are utilized by DESeq2's model?
  
    1. CPM (counts per million)  
    1. TPM (transcripts per kilobase million)  
    1. RPKM (reads per kilobase of exon per million reads/fragments mapped)  
    1. **median of ratios**
  
5. True or False: RPKM is the preferred normalization method to compare the expression of a given gene between samples.  
    1. True  
    1. **False**

6. True or False: The presence of differentially expressed genes between samples can alter the size factor (normalization factor).   
    1. True  
    1. **False**
  
7. Sample-level QC can be performed on normalized counts directly, but the log transformation of the normalized counts gives a better idea of similarity and differences
    1. **True**
    1. False
    
8. When I plot PC1 vs PC2 for my data, the samples are not separating along PC1 based on my experimental question/factor. Which of the following should I NOT do.
    1. Color my data points by other factors in the data systemtically to identify what is the major source of variation.
    1. I should plot other PCs to identify if the samples separate by the experimental factor along another PC
    1. I should remove some samples and recreate the PCA to see if it helps with the clustering
    1. **The data are too noisy, I need to throw out the experiment and start from scratch**

9. In the multiqc report one of my samples appeared to be fairly different from the others, it has different trends for the various plots. During sample-level QC, it is not clustering with the other replicates from the same sample group. What should I do?
    1. **At this point it is okay to remove it as an outlier. Once removed, redo the sample-level QC.**
    2. Don't remove it as an outlier
