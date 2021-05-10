#### Assignment 2 ####

#### Design formula
## Exercise
# 1. Suppose you wanted to study the expression differences between the two age groups in the metadata shown above, and major sources of variation were sex and treatment, how would the design formula be written?
# Ans: design = ~ sex + treatment + age

# 2. Based on our Mov10 metadata dataframe, which factors could we include in our design formula?
# Ans: design = ~ sampletype

# 3. What would you do if you wanted to include a factor in your design formula that is not in your metadata?
# Ans: Add that factor into your metadata.

## Exercise: How would the design formula be structured to perform the following analyses?
# 1. Test for the effect of treatment.
# Ans: design = ~ treatment

# 2. Test for the effect of genotype, while regressing out the variation due to treatment.
# Ans: design = ~ treatment + genotype

# 3. Test for the effect of genotype on the treatment effects.
# Ans: design = ~ genotype + treatment + genotype:treatment

#### Hypothesis testing
# 1. What is an appropriate hypothesis test if you are testing for expression differences across the developmental stages?
# Ans: Likelihood ratio test, because there are more than two groups for comparison.

# 2. Provide the line of code used to create the dds object.
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ sex + developmental_stage) # here we assume that the metadata includes two columns: sex and developmental_stage

# 3. Provide the line of code used to run DESeq2.
dds_lrt <- DESeq(dds, test="LRT", reduced = ~ sex) # since the full model is "sex + developmental_stage", the reduced model is then "sex". We don't need to use "1" here, because our reduced model has one factor.

# 4. The results of the differential expression analysis run identifies a group of genes that spike in expression between the first and second timepoints with no change in expression thereafter. How would we go about obtaining fold changes for these genes?
# Ans: We could use a Wald test to compare the groups we are interested in. It is not uncommon to run both tests, as you recall the LRT does not provide fold changes, and sometimes it can be helpful to further reduce down to a set of higher confidence genes by identifying those with higher fold change values.

#### Description of steps for DESeq2
# 1. Given the dispersion plot below, would you have any concerns regarding the fit of your data to the model?
# If not, what aspects of the plot makes you feel confident about your data?
# If so, what are your concerns? What would you do to address them?
# Ans: Yes, there are some concerns. The data does not scatter around the fitted curve, and the distribution of normalized counts are restricted in a small range. I would double check QC of my samples to make sure that there are no contamination or outliers.

#### Wald test results
# MOV10 Differential Expression Analysis: Control versus Knockdown
# Now that we have results for the overexpression results, do the same for the Control vs. Knockdown samples.

# 1. Create a contrast vector called contrast_kd.
contrast_kd <- c("sampletype", "MOV10_knockdown", "control")

# 2. Use contrast vector in the results() to extract a results table and store that to a variable called res_tableKD.
res_tableKD <- results(dds, contrast=contrast_kd, alpha = 0.05)

# 3. Shrink the LFC estimates using lfcShrink() and assign it back to res_tableKD.

res_tableKD <- lfcShrink(dds, coef="sampletype_MOV10_knockdown_vs_control", type="apeglm")

#### Summarizing results and extracting significant gene lists
# MOV10 Differential Expression Analysis: Control versus Knockdown

# 1. Using the same p-adjusted threshold as above (padj.cutoff < 0.05), subset res_tableKD to report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.
res_tableKD_tb <- res_tableKD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sigKD <- res_tableKD_tb %>%
  filter(padj < padj.cutoff)

# 2. How many genes are differentially expressed in the Knockdown compared to Control? How does this compare to the overexpression significant gene list (in terms of numbers)?
# Ans: There are 2,827 genes differentially expressed in the Knockdown compared to Control, and 4,774 genes differentially expressed in the Overexpression compared to Control. Therefore, less genes are present in the Knockdown significant gene list.

