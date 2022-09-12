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
  * [Count normalization](../lessons/02_DGE_count_normalization.md)
  * [Sample-level QC](../lessons/03_DGE_QC_analysis.md) (PCA and hierarchical clustering)

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
       <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Sometimes we are interested in how a gene changes over time. The Likelihood Ratio Test (LRT) is paricularly well-suited for this task.<br><br>This lesson will cover:<br>
             - Designing a LRT for a time-course analysis in DESeq2<br>
             - Identifying patterns that could be considered as diffentially expressed<br><br>
        </details>
    * [Gene annotation](../lessons/genomic_annotation.md)
        <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Next-generation analyses rely on annotations to provide a description for defining genes, transcripts and/or proteins. These annotations are often stored in publicly availible databases. <br><br>This lesson will cover:<br>
             - Describing these various annotation databaes<br>
             - Accessing annotations on these databases using R<br><br>
        </details>
    * [Functional analysis - over-representation analysis](../lessons/10_FA_over-representation_analysis.md)
        <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br>Oftentimes after completing an RNA-seq experiment, you will be left with a list of differentially expressed transcripts. You may be interested in knowing if these transcripts are enriched in certain pathways or find novel pathways. <br><br>This lesson will cover:<br>
             - Describing how functional enrichment tools yield statistically enriched functions or interactions<br>
             - Identify popular functional analysis tools for over-representation analysis<br><br>
        </details>
    * [Functional analysis - functional class scoring / GSEA](../lessons/11_FA_functional_class_scoring.md)
        <details>
       <summary><i>Click here for a preview of this lesson</i></summary>
         <br> While some functional analyses focus on large changes focused on a select few genes, functional class scoring (FCS) focuses on weaker but coordinated changes in sets of functionally related genes (i.e., pathways) that can also have significant effects. <br><br>This lesson will cover:<br>
             - Designing a GSEA analysis using GO and KEGG gene sets<br>
             - Evaluate the results of a GSEA analysis<br>
             - Discuss other tools and resources for identifying genes of novel pathways or networks<br><br>
        </details>

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
