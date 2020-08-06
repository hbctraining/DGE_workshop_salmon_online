# Workshop Schedule

## Pre-reading

1. [Workflow (raw data to counts)](../lessons/01a_RNAseq_processing_workflow.md)
1. [Experimental design considerations](https://hbctraining.github.io/Intro-to-rnaseq-hpc-salmon/lessons/experimental_planning_considerations.html)

## Day 1

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 10:30 | [Workshop Introduction](../lectures/Intro_to_workshop_dfhcc.pdf) | Radhika |
| 10:30 - 10:45 | [R refresher Q & A](../lessons/R_refresher.md) | Vishal |
| 10:45 - 11:15 | RNA-seq pre-reading discussion | Radhika/Vishal |
| 11:15 - 12:00 | [Intro to DGE / setting up DGE analysis](../lessons/01b_DGE_setup_and_overview.md) | Meeta |

### Self-Learning Part 1

1. [RNA-seq counts distribution](../lessons/01c_RNAseq_count_distribution.md)
1. [Count normalization](../lessons/02_DGE_count_normalization.md)
1. [Sample-level QC](../lessons/03_DGE_QC_analysis.md) (PCA and hierarchical clustering)

### Assignment #1
* All exercises from above lessons have been put together in [R script format](../homework/DGE_assignment_1.R).
* Add your solutions to the exercises in the downloaded .R file and upload the saved file to Dropbox on the Wednesday before the next class using [this link](https://www.dropbox.com/request/ud9JHqxu4OlfAT0X45k1).
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).
* **[Answer Key](../homework/DGE_assignment_1_answer_key.R)**

---

## Day 2

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:00 | Self-learning lessons discussion | Radhika |
| 11:00 - 11:30 | [Design formulas](../lessons/04a_design_formulas.md)  | Vishal |
| 11:30 - 12:00 | [Hypothesis testing and multiple test correction](../lessons/05a_hypothesis_testing.md) | Meeta |

### Self-Learning Part 2
1. [Description of steps for DESeq2](../lessons/04b_DGE_DESeq2_analysis.md)
1. [Wald test results](../lessons/05b_wald_test_results.md)
1. [Summarizing results and extracting significant gene lists](../lessons/05c_summarizing_results.md)
1. [Visualization](../lessons/06_DGE_visualizing_results.md)

### Assignment #2
* All exercises from above lessons have been put together in [R script format](../homework/DGE_assignment_2.R)
* Add your solutions to the exercises in the downloaded .R file and upload the saved file to Dropbox on the Monday before the next class using [this link](https://www.dropbox.com/request/c4Tt3OQC3nwCVEva3wBL).
* [Email us](mailto:hbctraining@hsph.harvard.edu) about questions regarding the homework that you need answered before the next class.
* Post questions that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945).

---

## Day 3

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:15 | Self-learning lessons discussion | Meeta/Vishal/Radhika |
| 11:15 - 12:00 | [Likelihood Ratio Test results](../lessons/08a_DGE_LRT_results.md) | Meeta |

### Self-Learning Part 3
1. [Time course analysis](../lessons/08b_time_course_analyses.md)
1. [Gene annotation](../lessons/genomic_annotation.md)
1. [Functional analysis - over-representation analysis](../lessons/10_FA_over-representation_analysis.md)
1. [Functional analysis - functional class scoring / GSEA](../lessons/11_FA_functional_class_scoring.md)

### Assignment #3
* There is no assignment submission for *Self-Learning Part 3*
* Post questions from the self learning materials that you would like to have reviewed in class [here](https://PollEv.com/hbctraining945) or [email us](mailto:hbctraining@hsph.harvard.edu).

---

## Day 4

| Time            |  Topic  | Instructor |
|:------------------------:|:------------------------------------------------:|:--------:|
| 10:00 - 11:00 | Questions about self-learning lessons | Radhika |
| 11:00 - 11:15 | [Summarizing workflow](../lessons/07_DGE_summarizing_workflow.md) | Vishal |
| 11:15 - 11:45 | Discussion, Q & A | All |
| 11:45 - 12:00 | [Wrap Up](../lectures/Workshop_wrapup.pdf) | Radhika |


## Resources
We have covered the inner workings of DESeq2 in a fair amount of detail such that when using this package you have a good understanding of what is going on under the hood. For more information on topics covered, we encourage you to take a look at the following resources:

* [DESeq2 vignette](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#theory-behind-deseq2)
* GitHub book on [RNA-seq gene level analysis](http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html)
* [Bioconductor support site](https://support.bioconductor.org/t/deseq2/) (posts tagged with `deseq2`) 

## Building on this workshop
* [Single-cell RNA-seq workshop](https://hbctraining.github.io/scRNA-seq/)
* [RMarkdown](https://hbctraining.github.io/Training-modules/Rmarkdown/)
* [Ggplot2 for functional analysis](https://hbctraining.github.io/Training-modules/Tidyverse_ggplot2/lessons/ggplot2.html)

## Other helpful links
* [Online learning resources](https://hbctraining.github.io/bioinformatics_online/lists/online_trainings.html)
* [All hbctraining materials](https://hbctraining.github.io/main)
