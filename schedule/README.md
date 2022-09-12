# Workshop Schedule

## Pre-reading

1. [Workflow (raw data to counts)](../lessons/01a_RNAseq_processing_workflow.md)
1. [Experimental design considerations](../lessons/experimental_planning_considerations.md)

## Day 1

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 10:30 | [Workshop Introduction](../lectures/Intro_to_workshop_all.pdf) | Meeta |
| 10:30 - 11:00 | RNA-seq pre-reading discussion | All |
| 11:00 - 11:45 | [Intro to DGE / setting up DGE analysis](../lessons/01b_DGE_setup_and_overview.md) | Jihe |
| 11:45 - 12:00 | Overview of self-learning materials and homework submission | All |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:
  * [RNA-seq counts distribution](../lessons/01c_RNAseq_count_distribution.md)
     <details>
          <summary><i>Click here for a preview of this lesson</i></summary>
            <br>Starting with the count matrix, we want to explore some characteristics of the RNA-seq data and evaluate the appropriate model to use. <br><br>This lesson will cover:<br>
                - Describing charactericstics of the RNA-seq count data<br>
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
   * **Copy over** your code from the exercises into a text file. 
   * **Upload the saved text file** to [Dropbox](https://www.dropbox.com/request/No9iPuFIJTL4oNJEpAZI) the **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

---

## Day 2

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:00 | Self-learning lessons discussion | All |
| 11:00 - 11:30 | [Design formulas](../lessons/04a_design_formulas.md)  | Jihe |
| 11:30 - 12:00 | [Hypothesis testing and multiple test correction](../lessons/05a_hypothesis_testing.md) | Meeta |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:
    * [Description of steps for DESeq2](../lessons/04b_DGE_DESeq2_analysis.md)
    * [Wald test results](../lessons/05b_wald_test_results.md)
    * [Summarizing results and extracting significant gene lists](../lessons/05c_summarizing_results.md)
    * [Visualization](../lessons/06_DGE_visualizing_results.md)

2. **Complete the exercises**:
   * Each lesson above contain exercises; please go through each of them.
   * **Copy over** your code from the exercises into a text file. 
   * **Upload the saved text file** to [Dropbox](https://www.dropbox.com/request/C7qCKdcXcyOBsB9m0NDj) the **day before the next class**.

### Questions?
* ***If you get stuck due to an error*** while runnning code in the lesson, [email us](mailto:hbctraining@hsph.harvard.edu) 
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

---

## Day 3

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:15 | Self-learning lessons discussion | Meeta |
| 11:15 - 12:00 | [Likelihood Ratio Test results](../lessons/08a_DGE_LRT_results.md) | Jihe |

### Before the next class:

1. Please **study the contents** and **work through all the code** within the following lessons:
    * [Time course analysis](../lessons/08b_time_course_analyses.md)
    * [Gene annotation](../lessons/genomic_annotation.md)
    * [Functional analysis - over-representation analysis](../lessons/10_FA_over-representation_analysis.md)
    * [Functional analysis - functional class scoring / GSEA](../lessons/11_FA_functional_class_scoring.md)

2. **There is no assignment submission**

### Questions?
* Post any **conceptual questions** that you would like to have **reviewed in class** [here](https://PollEv.com/hbctraining945).

---

## Day 4

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:00 | Questions about self-learning lessons | All |
| 11:00 - 11:15 | [Summarizing workflow](../lessons/07_DGE_summarizing_workflow.md) | Jihe |
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
* [Functional analysis visualization](https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html)

## Building on this workshop
* [Single-cell RNA-seq workshop](https://hbctraining.github.io/scRNA-seq/)
* [RMarkdown](https://hbctraining.github.io/Training-modules/Rmarkdown/)
* [ggplot2 for functional analysis](https://hbctraining.github.io/Training-modules/Tidyverse_ggplot2/lessons/03_ggplot2.html)

## Other helpful links
* [Online learning resources](https://hbctraining.github.io/bioinformatics_online/lists/online_trainings.html)
* [All hbctraining materials](https://hbctraining.github.io/main)
