---
title: "From raw sequence reads to count matrix"
author: "Meeta Mistry, Radhika Khetani, Mary Piper"
date: "May 12, 2017"
---

Approximate time: 

## Learning Objectives 

* Detailed description of dataset
* Numbered descriptions of workflow steps leading up to count matrix
* Link out to Illumina video


## RNA-seq processing workflow

Before we get started with differential gene expression, it's important to know where the count matrix came from. In this lesson we will briefly discuss the RNA-processing pipeline and the **different steps we take to go from raw sequencing reads to count matrix**. 

<p align="center">
<img src="../img/workflow-salmon-DGE.png" >
</p>


### 1. RNA Extraction and library preparation

Before RNA can be sequenced, it must first be extracted and separated from its cellular environment prepared into a cDNA library. There are a number of steps involved which are outlined in the figure below, and in parallel there are various quality checks implemented to make sure wehave quality RNA to move foward with. We briefly describe some of these steps, but also encourage you to access the resources linked at the end of this lesson for more detailed information.

**a. Enriching for RNA.** Once the sample has been treated with DNAse to remove any contaminating DNA sequence, the sample undergoes either selection of the mRNA (polyA selection) or depletion of the rRNA. 

> **RNA Quality check**: Traditionally, RNA integrity was assessed via gel electrophoresis by visual inspection of the ribosomal RNA bands. Since it is subjective and has been shown to be inconsistent, Agilent developed a software algorithm that allows for the calculation of an RNA Integrity Number (RIN) which facilitates the interpretation and reproducibility, of RNA quality assessments. RIN provides a means by which samples can be compared in a standardized manner.

**b. **



<p align="center">
<img src="../img/library_prep.png" >
</p>

*Image source: [Introduction to differential gene expression analysis using RNA-seq](https://chagall.med.cornell.edu/RNASEQcourse/Intro2RNAseq.pdf)*

### 2. Seqeuncing (Illumina)

### 3. Quality control of raw sequencing data

### 4. Mapping reads and quantification

### 5. Quality control of mapped sequence reads
