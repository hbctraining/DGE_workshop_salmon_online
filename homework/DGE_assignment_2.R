#### Assignment 2 ####

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

