#### Assignment 2 ####

#### Design formula
# 1. How would the design formula be structured to perform the following analyses?

# 2. Test for the effect of treatment.

# 3. Test for the effect of genotype, while regressing out the variation due to treatment.

# 4. Test for the effect of genotype on the treatment effects.

#### Hypothesis testing
# 1. What is an appropriate hypothesis test if you are testing for expression differences across the developmental stages?

# 2. Provide the line of code used to create the dds object.

# 3. Provide the line of code used to run DESeq2.

# 4. The results of the differential expression analysis run identifies a group of genes that spike in expression between the first and second timepoints with no change in expression thereafter. How would we go about obtaining fold changes for these genes?
  
#### Description of steps for DESeq2
# 1. Given the dispersion plot below, would you have any concerns regarding the fit of your data to the model?
# If not, what aspects of the plot makes you feel confident about your data?
# If so, what are your concerns? What would you do to address them?

#### Wald test results
# MOV10 Differential Expression Analysis: Control versus Knockdown
# Now that we have results for the overexpression results, do the same for the Control vs. Knockdown samples.

# 1. Create a contrast vector called contrast_kd.

# 2. Use contrast vector in the results() to extract a results table and store that to a variable called res_tableKD.

# 3. Shrink the LFC estimates using lfcShrink() and assign it back to res_tableKD.

#### Summarizing results and extracting significant gene lists
# MOV10 Differential Expression Analysis: Control versus Knockdown

# 1. Using the same p-adjusted threshold as above (padj.cutoff < 0.05), subset res_tableKD to report the number of genes that are up- and down-regulated in Mov10_knockdown compared to control.

# 2. How many genes are differentially expressed in the Knockdown compared to Control? How does this compare to the overexpression significant gene list (in terms of numbers)?

