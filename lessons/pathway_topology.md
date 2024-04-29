### SPIA
The [SPIA (Signaling Pathway Impact Analysis)](http://bioconductor.org/packages/release/bioc/html/SPIA.html) tool can be used to integrate the lists of differentially expressed genes, their fold changes, and pathway topology to identify affected pathways. The blog post from [Getting Genetics Done](https://gettinggeneticsdone.blogspot.com/2012/03/pathway-analysis-for-high-throughput.html) provides a step-by-step procedure for using and understanding SPIA. 

```r
# Set-up

BiocManager::install("SPIA")
library(SPIA)

## Significant genes is a vector of fold changes where the names are ENTREZ gene IDs. The background set is a vector of all the genes represented on the platform.

background_entrez <- res_entrez$entrezid

sig_res_entrez <- res_entrez[which(res_entrez$padj < 0.05), ]

sig_entrez <- sig_res_entrez$log2FoldChange

names(sig_entrez) <- sig_res_entrez$entrezid

head(sig_entrez)
```


Now that we have our background and significant genes in the appropriate format, we can run SPIA:

```r

spia_result <- spia(de=sig_entrez, all=background_entrez, organism="hsa")

head(spia_result, n=20)
```

SPIA outputs a table showing significantly dysregulated pathways based on over-representation and signaling perturbations accumulation. The table shows the following information: 

- `pSize`: the number of genes on the pathway
- `NDE`: the number of DE genes per pathway
- `tA`: the observed total perturbation accumulation in the pathway
- `pNDE`: the probability to observe at least NDE genes on the pathway using a hypergeometric model (similar to ORA)
- `pPERT`: the probability to observe a total accumulation more extreme than tA only by chance
- `pG`: the p-value obtained by combining pNDE and pPERT
- `pGFdr` and `pGFWER` are the False Discovery Rate and Bonferroni adjusted global p-values, respectively
- `Status`: gives the direction in which the pathway is perturbed (activated or inhibited)
- `KEGGLINK` gives a web link to the KEGG website that **displays the pathway image** with the differentially expressed genes highlighted in red

We can view the significantly dysregulated pathways by viewing the over-representation and perturbations for each pathway.

```r
plotP(spia_result, threshold=0.05)
```

![perturbed_pathway](../img/spia_plot.png)

In this plot, each pathway is a point and the coordinates are the log of pNDE (using a hypergeometric model) and the p-value from perturbations, pPERT. The oblique lines in the plot show the significance regions based on the combined evidence.

If we choose to explore the significant genes from our dataset occurring in these pathways, we can subset our SPIA results:

```r
## Look at pathway 03013 and view kegglink
subset(spia_result, ID == "03013")
```

Then, click on the KEGGLINK, we can view the genes within our dataset from these perturbed pathways:

![perturbed_pathway](../img/hsa05222.png)
