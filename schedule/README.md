# Workshop Schedule

## Pre-reading

1. [Workflow (raw data to counts)](../lessons/01a_RNAseq_processing_workflow.md)
1. [Experimental design considerations](../lessons/experimental_planning_considerations.md)

## Day 1

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 10:30 | [Workshop Introduction](../lectures/Intro_to_workshop_all.pdf) | Meeta |
| 10:30 - 11:00 | RNA-seq pre-reading discussion | All |
| 11:00 - 11:45 | [Intro to DGE / setting up DGE analysis](../lessons/01b_DGE_setup_and_overview.md) | Noor |
| 11:45 - 12:00 | Overview of self-learning materials and homework submission | Meeta |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:

  * [RNA-seq counts distribution](../lessons/01c_RNAseq_count_distribution.md)
     <details>
          <summary><i>Click here for a preview of this lesson</i></summary>
            <br>Starting with the count matrix, we want to explore some characteristics of the RNA-seq data and evaluate the appropriate model to use. <br><br>This lesson will cover:<br>
                - Describing characteristics of the RNA-seq count data<br>
                - Understanding different statistical methods to model the count data<br>
                - Explaining the benefits of biological replicates<br><br>
           </details>
        
  * [Count normalization](../lessons/02_DGE_count_normalization.md)
     <details>
          <summary><i>Click here for a preview of this lesson</i></summary>
            <br>Count normalization is an import data pre-processing step before the differential expression analysis. <br><br>This lesson will cover:<br>
                - Describing "uninteresting factors" to consider during normalization<br>
                - Understanding different normalization methods and their corresponding use cases<br>
                - Generating a matrix of normalized counts using DESeq2's median of ratios method<br><br>
           </details>
  
  * [Sample-level QC](../lessons/03_DGE_QC_analysis.md) (PCA and hierarchical clustering)
     <details>
          <summary><i>Click here for a preview of this lesson</i></summary>
            <br>Next, we want to check the quality of count data, to make sure that the samples are good. 
            <br><br>This lesson will cover:<br>
                - Understanding the importance of similarity analysis between samples<br>
                - Describing Principal Component Analysis (PCA) and interpreting PCA plots from RNA-seq data<br>
                - Performing hierarchical clustering and plotting correlation metrics<br><br>
           </details>

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your solutions into the [Google Form](https://docs.google.com/forms/d/e/1FAIpQLScUcYzyiM_dAsgQdNx9ECzCX3lKrTHTwmUKux9u8VyP2JDLNQ/viewform?usp=sf_link) **the day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

---

## Day 2

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:00 | Self-learning lessons discussion | All |
| 11:00 - 11:30 | [Design formulas](../lessons/04a_design_formulas.md)  | Noor |
| 11:30 - 12:00 | [Hypothesis testing and multiple test correction](../lessons/05a_hypothesis_testing.md) | Meeta |

### Before the next class:

I. Please **study the contents** and **work through all the code** within the following lessons:
   1. [Description of steps for DESeq2](../lessons/04b_DGE_DESeq2_analysis.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br> The R code required to perform differential gene expression analysis is actually quite simple. Running the `DESeq()` function will carry out the various steps involved. It is important that you have some knowledge of what is happening under the hood, to be able to fully understand and interpret the results <br><br>In this lesson you will:<br>
             - Examine size factors and learn about sources that cause observed variation in values <br>
             - Explore the gene-wise dispersion estimates as they relate back the mean-variance relationship <br>
             - Critically evaluate a dispersion plot <br><br>
        </details>

   2. [Wald test results](../lessons/05b_wald_test_results.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br> We have run the analysis, and now it's time to explore the results!  <br><br>In this lesson you will:<br>
             - Learn how to extract results for specific group comparisons <br>
             - Explore the information presented in the results table (different statistics and their importance) <br>
             - Understand the different levels of filtering that are applied in DESeq2 by default (and why they are important) <br><br>
        </details>
        
        
   3. [Summarizing results and extracting significant gene lists](../lessons/05c_summarizing_results.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br> Once you have your results, it is useful to summarize the information. Here, we get a snapshot of the number of differentially expressed genes that are identified from the different comparisons. <br><br>
        </details>
        
 4. [Visualization](../lessons/06_DGE_visualizing_results.md)
      <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>A picture is worth a thousand words. In our case, a figure is worth a thousand (or 30 thousand) data points. When working with large scale data, it can be helpful to visualize results and get a big picture perspective of your findings. <br><br>In this lesson you will:<br>
            - Explore different plots for data visualization <br>
            - Create a volcano plot to evaluate the relationship between different statistics from the results table <br>
            - Create a heatmap for visualization of differentially expressed genes <br><br>
        </details>

II. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your solutions into the [Google Form](https://docs.google.com/forms/d/e/1FAIpQLSfVELkIcVN4wyJ2aNrowgxiuat5uUXCXACj8QN4MfTK5Yr-Zw/viewform?usp=sf_link) **the day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

---

## Day 3

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:15 | Self-learning lessons discussion | All |
| 11:15 - 12:00 | [Likelihood Ratio Test results](../lessons/08a_DGE_LRT_results.md) | Meeta |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:
    * [Time course analysis](../lessons/08b_time_course_analyses.md)
       <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Sometimes we are interested in how a gene changes over time. The Likelihood Ratio Test (LRT) is paricularly well-suited for this task.<br><br>This lesson will cover:<br>
             - Designing a LRT for a time-course analysis in DESeq2<br>
             - Identifying patterns in our list of differentially expressed genes<br><br>
        </details>
    * [Gene annotation](../lessons/genomic_annotation.md)
        <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Next-generation analyses rely on annotations to provide a description for defining genes, transcripts and/or proteins. These annotations are often stored in publicly available databases. <br><br>This lesson will cover:<br>
             - Describing the various annotation databases<br>
             - Accessing annotations from one of these databases using R<br><br>
        </details>
    * [Functional analysis - over-representation analysis](../lessons/10_FA_over-representation_analysis.md)
        <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Oftentimes after completing an RNA-seq experiment, you will be left with a list of differentially expressed transcripts. You may be interested in knowing if these transcripts are enriched in certain biologically-relevant contexts. <br><br>This lesson will cover:<br>
             - Describing how functional enrichment tools yield statistically enriched functional categories or interactions<br>
             - Identifying enriched Gene Ontology terms using the R package, clusterProfiler <br><br>
        </details>
    * [Functional analysis - functional class scoring / GSEA](../lessons/11_FA_functional_class_scoring.md)
        <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br> While some functional analyses focus on large changes focused on a select few genes, functional class scoring (FCS) focuses on weaker but coordinated changes in sets of functionally related genes (i.e., pathways) that can also have significant effects. <br><br>This lesson will cover:<br>
             - Designing a GSEA analysis using GO and/or KEGG gene sets<br>
             - Evaluating the results of a GSEA analysis<br>
             - Discussing other tools and resources for identifying genes of novel pathways or networks<br><br>
        </details>

2. **There is no assignment submission**

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 

---

## Day 4

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:00 | Questions about self-learning lessons | All |
| 11:00 - 11:15 | [Summarizing workflow](../lessons/07_DGE_summarizing_workflow.md) | Noor |
| 11:15 - 11:45 | Discussion, Q & A | All |
| 11:45 - 12:00 | [Wrap Up](../lectures/Workshop_wrapup_all.pdf) | Meeta |

## Answer keys
* [Day 1 Answer Key](../homework/DGE_assignment_1_answer_key.R)
* [Day 2 Answer Key](../homework/DGE_assignment_2_answer_key.R)

## Resources
We have covered the inner workings of DESeq2 in a fair amount of detail such that when using this package you have a good understanding of what is going on under the hood. For more information on topics covered, we encourage you to take a look at the following resources:

* [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory-behind-deseq2)
* GitHub book on [RNA-seq gene level analysis](http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html)
* [Bioconductor support site](https://support.bioconductor.org/t/deseq2/) (posts tagged with `deseq2`) 
* [Enrichment analysis book](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html)
   * [Visualization: Functional (Enrichment) analysis](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)

## Building on this workshop
* [Single-cell RNA-seq workshop](https://hbctraining.github.io/scRNA-seq/)
* [RMarkdown](https://hbctraining.github.io/Training-modules/Rmarkdown/)
* [ggplot2 for functional analysis](https://hbctraining.github.io/Training-modules/Tidyverse_ggplot2/lessons/03_ggplot2.html)

## Other helpful links
* [Online learning resources](https://hbctraining.github.io/bioinformatics_online/lists/online_trainings.html)
* [All hbctraining materials](https://hbctraining.github.io/main)
