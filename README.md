# TADA-A
A statistical framework for mapping risk genes from *de novo* mutations in whole-genome sequencing studies

## 1. Introduction
With the fast pace of technology revolution in the field of genomics, whole-exome sequencing (WES) and whole-genome sequencing (WGS) have become more and more affordable. We have seen rapid accumulation of these sequencing data from cohorts with different traits and diseases and are about to see even more in the near future. Now, the major challenge is what we could learn from these data. Specifically, **could we use mutation data to learn what types of mutations are relevant to disease eitiology and which genes are likely to be risk genes**. TADA-A is a statistical framework designed to answer these important questions using *de novo mutations* (DNM), one type of mutations that spontaneously arise in an offspring and not present in its parents. The two main types of input data for TADA-A is DNM and functional/conservational annotation data. It works by integrating DNM information across different studies while accounting for technical variations that might affect observed DNM mutation rates among these studies. It first adjusts baseline mutation rates for each study and estimates the effect sizes of annotations using DNM data from all studies. Finally, it uses relevant annotations to predict the deleteriouss of each mutation and predict disease risk genes. 

## 2. Prerequisites

### 2.1 bedtools
Bedtools need to be installed and added to your PATH. We suggest using v2.17.0, which has been provided in the companion files for you to download. Other versions might be incompatible because of modification of input arguments for some of the sub-functions of bedtools. 

### 2.2 bigWigAverageOverBed
This executable has been added in the `external_tools` folder.

## 3. 
