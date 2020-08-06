#### Assignment 1 ####

#### RNA-seq counts distribution
# 1. Evaluate the relationship between mean and variance for the control replicates (Irrel_kd samples). Note the differences or similarities in the plot compared to the one using the overexpression replicates.

mean_counts_ctrl <- apply(data[,1:3], 1, mean)        #select column 1 to 3, which correspond to Irrel_kd samples
variance_counts_ctrl <- apply(data[,1:3], 1, var)
df_ctrl <- data.frame(mean_counts_ctrl, variance_counts_ctrl)
ggplot(df_ctrl) +
  geom_point(aes(x=mean_counts_ctrl, y=variance_counts_ctrl)) + 
  scale_y_log10(limits = c(1,1e9)) +
  scale_x_log10(limits = c(1,1e9)) +
  geom_abline(intercept = 0, slope = 1, color="red")

# Ans: The plot of mean and variance for the control replicates is similar to that of overexpression replicates shown in the lesson.

# 2. An RNA-seq experiment was conducted on mice forebrain to evaluate the effect of increasing concentrations of a treatment. For each of the five different concentrations we have n = 5 mice for a total of 25 samples. If we observed little to no variability between replicates, what might this suggest about our samples?

# Ans: The lack of variability between replicates suggests that we are possibly dealing with technical replicates. With true biologcal replicates we expect some amount of variability. If you have technical replicates, you do not want to be using DESeq2 because we will be using the NB to account for overdispersion, which doesn't exist.

# 3. What type of mean-variance relationship would you expect to see for this dataset?
# Ans: mean == variance. A Poisson would be more appropriate.

#### Count normalization
# 1. Suppose we have sample names matching in the counts matrix and metadata file, but they are in different order. Write the line(s) of code to create a new matrix with columns re-ordered such that they are identical to the row names of the metadata.

idx <- match(rownames(meta), colnames(data))
data_reordered <- data[, idx]

## OR

data_reordered <- data[, match(rownames(meta), colnames(data))]

## OR

txi$counts <- txi$counts[, match(rownames(meta), colnames(txi$counts))]

#### Sample-level QC
# 1. What does the above plot tell you about the similarity of samples?
# Ans: Samples from different experimental groups are different, while replicates within the same group are similar.

# 2. Does it fit the expectation from the experimental design?
# Ans: Yes, it does.

# 3. What do you think the %variance information (in the axes titles) tell you about the data in the context of the PCA?
# Ans: PC1 is associated with 37% of the variance in the data, and PC2 is associated with 16%. 
# *You won't see the %variance information when you make PCA plots without the plotPCA() function. But, there are tools that let you explore this.*
