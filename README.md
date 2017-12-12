# TADA-A
A statistical framework for mapping risk genes from *de novo* mutations in whole-genome sequencing studies

## 1. Introduction
With the fast pace of technology revolution in the field of genomics, whole-exome sequencing (WES) and whole-genome sequencing (WGS) have become more and more affordable. We have seen rapid accumulation of these sequencing data from cohorts with different traits and diseases and are about to see even more in the near future. Now, the major challene is what we could learn from these data. Specifically, **could we use mutation data to learn what types of mutations are relevant to disease eitiology and which genes are likely to be risk genes**. TADA-A is a statistical framework designed to answer these important questions using *de novo mutations* (DNM), one type of mutations that spontaneously arise in an offspring and not present in its parents. 

## 2. Prerequisites

### 2.1 bedtools
Bedtools need to be installed and added to your PATH. We suggest using v2.17.0, which has been provided in the companion files for you to download. Other versions might be incompatible because of modification of input arguments for some of the sub-functions of bedtools. 

### bigWigAverageOverBed
This executable has been added in the `external_tools` folder.


