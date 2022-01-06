---
  title: "Functional Analysis of Differential Expression Results "
author: "Analyst name"
contact: analyst_email@hsph.harvard.edu
project: PI Bulk RNA-seq Analysis
date: "`r Sys.Date()`"
output:
  html_document:
  code_folding: hide
df_print: paged
highlights: pygments
number_sections: true
self_contained: true
theme: cosmo
toc: true
toc_depth: 2
toc_float:
  collapsed: true
smooth_scroll: true
fig_width: 6
fig_height: 5
---
  
  # Overview
  
  - **Principal Investigator:** PI name
- **Researcher:** Researcher name - often the person of main contact
- **Experiment:** Experiment description in one sentence
- **Experimental details:** Important details of the experiment, often acquired during the initial meeting. Could have subsections, such as:
  
  - **Hypotheses:**
  - **Experimental Design:**
  - **Experimental Goals:**
  - **Expectations:**
  
  -**Functional analysis:** This report is specifically exploring ...

```{r setup, echo = FALSE, cache = FALSE}
knitr::opts_chunk$set(dev = c('png', 'cairo_pdf'),
                      fig.align = 'center', 
                      fig.height = 5, 
                      fig.width = 7,
                      pdf.options(encoding = "ISOLatin9.enc"),
                      fig.path='figures/',
                      warning=FALSE, 
                      message=FALSE,
                      cache = FALSE,
                      dev = c("png", "pdf"),
                      error = TRUE,
                      highlight = TRUE,
                      prompt = FALSE,
                      tidy = FALSE)
```

```{r setup, message=FALSE, warning=FALSE}
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(cowplot)
library(png)
library(grid)
library(pheatmap)
library(RColorBrewer)
library(magick)
library(AnnotationHub)
library(ensembldb)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DEGreport)
library(png)

# Set ggplot2 default theme
ggplot2::theme_set(theme_light(base_size = 14))
```

# Functional analysis {.tabset}

We have performed functional analysis using over-representation analysis (ORA) and gene set enrichment analysis (GSEA) methods. We have performed each of these methods using gene ontology (GO) and Kyoto Encyclopedia of Genes and Genomes (KEGG) gene sets.

```{r fxl_setup, message=FALSE, warning=FALSE}
# Convert Ensembl gene identifiers to Entrez
library(AnnotationHub)
library(ensembldb)

allGenes <- ann_results$Row.names %>%
  as.character()
sigGenes <- as.character(sig_results$Row.names)


# Get gene annotations based on reference data version used for alignment/quantification from BiomaRt

# Specify the Ensembl release ## dataset used by bcbio (check archives if needed)

# List current versions of datasets
#ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL")
#datasets <- listDatasets(ensembl)

# Identify the proper archive to use for corresponding Ensembl release
#archives <- listEnsemblArchives()

# This is example code for using the Ensembl 99 release for the human genome
ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                   dataset = "hsapiens_gene_ensembl",
                   host = "jan2020.archive.ensembl.org")

## Build biomaRt query
# filters = listFilters(ensembl)
# attributes = listAttributes(ensembl)

gene_annotations <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_id', 'gene_biotype', 'external_gene_name', 'description'), 
                          filters = 'ensembl_gene_id', 
                          values = rownames(allGenes), 
                          mart = ensembl)

```

```{r sig_tables, message=FALSE, warning=FALSE}
sigResults <- as.data.frame(ann_results)[ann_results$Row.names %in% sigGenes, ]

foldChanges <- sigResults$log2FoldChange
names(foldChanges) <- sigResults$Row.names

fdrOrdered <- ann_results %>%
  .[order(.$padj), ]

```


## GO over-representation analysis

Gene Ontology (GO) term over-representation analysis is a technique for interpreting sets of genes making use of the Gene Ontology system of classification, in which genes are assigned to a set of predefined bins depending on their functional characteristics.

There were ## over-represented GO terms identified.
Many of the most significant terms are related to ..., among others.

```{r enrich_go, message=FALSE, warning=FALSE}
library(clusterProfiler)
library(org.Mm.eg.db)
# Run GO enrichment analysis
ego <- enrichGO(
  sigGenes,
  "ENSEMBL",
  universe = allGenes,
  OrgDb = org.Mm.eg.db,
  ont = "ALL",
  pAdjustMethod = "BH",
  qvalueCutoff = 0.05,
  readable = TRUE)
saveRDS(ego, file = "data/")

# Show summary data.frame
egoSummary <- ego@result
write.csv(egoSummary,
          file = "results/", quote = FALSE)
knitr::kable(egoSummary[1:20, ]) 
```

[Download GO ORA results](./results/)

### DOTplot

```{r dotplot, fig.width=8, fig.height=12, message=FALSE, warning=FALSE}
# Dotplot of top 25
dotplot(ego, showCategory = 25)
```

### Category netplot

```{r cnetplot, fig.width=15, fig.height=12, message=FALSE, warning=FALSE}
# Cnet plot with genes colored by fold changes for top 5 most significant GO
# processes
foldChanges <- ifelse(foldChanges < -2, -2, foldChanges)
foldChanges <- ifelse(foldChanges > 2, 2, foldChanges)

ego2 <- ego
ego2@result <- ego@result[c(), ]

cnetplot(ego2,
         categorySize = "pvalue",
         showCategory = 8,
         foldChange = foldChanges,
         vertex.label.cex = 0.5)

```



## KEGG over-representation analysis

KEGG pathways were also tested for over-representation, and ## pathways were identified, including ... We used a p-adjusted value cut-off of 0.05 significance level. The significantly enriched pathways are displayed below, and only the differentially expressed genes are colored (the intensity of shading does not have any significance).

```{r enrich_kegg, message=FALSE, warning=FALSE}
library(pathview)

# Extract the Entrez IDs
allEntrez <- allGenes_ids$entrezid
allEntrez <- unlist(allGenes_ids$entrezid)

# Remove all NA values and only keep unique IDs
allEntrez <- as.character(unique(allEntrez[which(!(is.na(allEntrez)))]))
allEntrez_df <- allGenes_ids[allGenes_ids$entrezid %in% allEntrez, ]

# Combine the results for all genes with the Entrez ID
allEntrez_res <- merge(allEntrez_df, data.frame(ann_results), by.x="gene_id", by.y="Row.names")

allEntrez_res <- allEntrez_res[which(!(duplicated(allEntrez_res$entrezid))), ]

# Named list of log2 foldchanges
allEntrez_fc <- allEntrez_res$log2FoldChange
names(allEntrez_fc) <- allEntrez_res$entrezid

# Sort by fold change values
allEntrez_fc <- sort(allEntrez_fc, decreasing = TRUE)
sigEntrez <- unlist(allGenes_ids[allGenes_ids$gene_id %in% sigGenes, "entrezid"])
sigEntrez <- as.character(unique(sigEntrez[which(!(is.na(sigEntrez)))]))


kegg <- enrichKEGG(
  gene = sigEntrez,
  universe = allEntrez,
  organism = "mmu")
saveRDS(kegg, file = "data/")

# Show KEGG summary data.frame
keggSummary <- kegg@result

sigKegg <- keggSummary[keggSummary$p.adjust < 0.05, ]

write.csv(sigKegg, 
          file = "results/")

knitr::kable(sigKegg)

pathways <- sigKegg$ID

sigEntrez_fc <- allEntrez_fc[which(names(allEntrez_fc) %in% sigEntrez)]

# setwd("results/kegg_pathways/ora/")
# 
# for (pathway in pathways){
# pathview(gene.data = sigEntrez_fc,
#             pathway.id = pathway,
#             species = "mmu",
#             limit = list(gene = 2, cpd = 1))
# }
# 
# setwd("../../")

figures <- list.files("results/kegg_pathways/ora/", pattern = "pathview", full.names = TRUE)

figure_idx <- map_int(seq_along(pathways), function(a) {
  which(str_detect(pattern = pathways[a], string = figures))
})

img <- map(seq_along(figures), function(a){
  image_read(list.files("results/kegg_pathways/ora", pattern = "pathview", full.names = TRUE)[a]) %>%
    image_resize("570x380") %>%
    image_colorize(35, "white")
})

map(figure_idx, function(a){
  ggdraw() +
    draw_image(img[[a]])
})

keggPlotsDir <- "results/kegg_pathways/ora/"
for (file in list.files(keggPlotsDir, pattern = "pathview")){
  img <- image_read(paste0(keggPlotsDir, "/", file))
  print(ggdraw() +
          draw_image(img))
}
```

[Download KEGG ORA results](./results/)

## GO GSEA analysis

A common approach in analyzing gene expression profiles was identifying differential expressed genes that are deemed interesting. The enrichment analysis we demonstrated previously were based on these differentially expressed genes. This approach will find genes where the difference is large, but it will not detect a situation where the difference is small. Gene Set Enrichment Analysis (GSEA) directly addresses this limitation. All genes can be used in GSEA; GSEA aggregates the per gene statistics across genes within a gene set, therefore making it possible to detect situations where all genes in a predefined set change in a small but coordinated way. Since it is likely that many relevant phenotypic differences are manifested by small but consistent changes in a set of genes.

We are using the log2 fold changes as input. By using the log2 fold changes as the input, we are identifying processes with genes that exhibit coordinated fold changes that are larger than might be expected by chance. 

```{r gse_go, message=FALSE, warning=FALSE}
# Now run GSEA
gse <- gseGO(
  geneList = allEntrez_fc,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "ENTREZID",
  nPerm = 1000,
  minGSSize = 100,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = FALSE)
saveRDS(gse, file = "data/")

# Write out summary data.frame
gseSummary <- gse@result

write.csv(gseSummary,
          file = "results/",
          quote = FALSE)
knitr::kable(gseSummary[1:30, ])
```

[Download GO GSEA results](./results/)


## KEGG GSEA analysis

We can also perform GSEA analysis with clusterProfiler using KEGG gene sets. With this analysis we identified the ## enriched pathways shown in the table below. In these images we are showing the foldchanges for all genes tested for differential expression.

```{r kegg_gsea, message=FALSE, warning=FALSE}

# GSEA using gene sets from KEGG pathways
gseaKEGG <- gseKEGG(
  geneList = allEntrez_fc,
  organism = "mmu",
  nPerm = 2000,
  minGSSize = 10,
  pvalueCutoff = 0.05,
  verbose = FALSE)
saveRDS(gseKEGG, file = "data/")

# Extract the GSEA results
gseaKEGGSummary <- slot(gseaKEGG, "result") %>% as_tibble()
write.csv(gseaKEGGSummary,
          file = "results/",
          quote = FALSE)

knitr::kable(gseaKEGGSummary)
# dplyr must be unloaded at this step for pathview to work
suppressWarnings(detach("package:dplyr", unload = TRUE))

# If there is an error at this step, there may be a pathway that is not found by
# pathview package. In this case, you may need to run the pathview command above
# by specifying the index of the pathways you would like to print out in place
# of `x`.
pathways <- gseaKEGGSummary$ID

current <- getwd()

setwd("results/kegg_pathways/gsea")

for (pathway in pathways){
  pathview(gene.data = allEntrez_fc,
           pathway.id = pathway,
           species = "mmu",
           limit = list(gene = 2, cpd = 1))
}


setwd(current)

figures <- list.files("results/kegg_pathways/gsea", pattern = "pathview", full.names = TRUE)

figure_idx <- map_int(seq_along(pathways), function(a) {
  which(str_detect(pattern = pathways[a], string = figures))
})

img <- map(seq_along(figures), function(a){
  image_read(list.files("results/kegg_pathways/gsea", pattern = "pathview", full.names = TRUE)[a]) %>%
    image_resize("570x380") %>%
    image_colorize(35, "white")
})

map(figure_idx, function(a){
  ggdraw() +
    draw_image(img[[a]])
})

keggPlotsDir <- "results/kegg_pathways/gsea"
for (file in list.files(keggPlotsDir, pattern = "pathview")){
  img <- image_read(paste0(keggPlotsDir, "/", file))
  print(ggdraw() +
          draw_image(img))
}

```

[Download KEGG GSEA results](./results/)

# Session information

```{r sessionInfo}
sessionInfo()
```