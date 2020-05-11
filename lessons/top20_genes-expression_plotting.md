#### Using `ggplot2` to plot multiple genes (e.g. top 20)

Often it is helpful to check the expression of multiple genes of interest at the same time. This often first requires some data wrangling.

We are going to plot the normalized count values for the **top 20 differentially expressed genes (by padj values)**. 

To do this, we first need to determine the gene names of our top 20 genes by ordering our results and extracting the top 20 genes (by padj values):

```r
## Order results by padj values
top20_sigOE_genes <- res_tableOE_tb %>% 
        arrange(padj) %>% 	#Arrange rows by padj values
        pull(gene) %>% 		#Extract character vector of ordered genes
        head(n=20)		#Extract the first 20 genes
```

Then, we can extract the normalized count values for these top 20 genes:

```r
## normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
        filter(gene %in% top20_sigOE_genes)
```

Now that we have the normalized counts for each of the top 20 genes for all 8 samples, to plot using `ggplot()`, we need to gather the counts for all samples into a single column to allow us to give ggplot the one column with the values we want it to plot.

The `gather()` function in the **tidyr** package will perform this operation and will output the normalized counts for all genes for *Mov10_oe_1* listed in the first 20 rows, followed by the normalized counts for *Mov10_oe_2* in the next 20 rows, so on and so forth.

<img src="../img/melt_wide_to_long_format.png" width="800">

```r
# Gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
  gather(colnames(top20_sigOE_norm)[2:9], key = "samplename", value = "normalized_counts")

## check the column header in the "gathered" data frame
View(gathered_top20_sigOE)
```

Now, if we want our counts colored by sample group, then we need to combine the metadata information with the melted normalized counts data into the same data frame for input to `ggplot()`:

```r
gathered_top20_sigOE <- inner_join(mov10_meta, gathered_top20_sigOE)
```

The `inner_join()` will merge 2 data frames with respect to the "samplename" column, i.e. a column with the same column name in both data frames.

Now that we have a data frame in a format that can be utilised by ggplot easily, let's plot! 

```r
## plot using ggplot2
ggplot(gathered_top20_sigOE) +
        geom_point(aes(x = symbol, y = normalized_counts, color = sampletype)) +
        scale_y_log10() +
        xlab("Genes") +
        ylab("log10 Normalized Counts") +
        ggtitle("Top 20 Significant DE Genes") +
        theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	theme(plot.title = element_text(hjust = 0.5))
```

<img src="../img/sig_genes_gather_salmon.png" width="700">
