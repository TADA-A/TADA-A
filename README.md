# TADA-A
A statistical framework for mapping risk genes from *de novo* mutations in whole-genome sequencing studies

## 1. Introduction
With the fast pace of technology revolution in the field of genomics, whole-exome sequencing (WES) and whole-genome sequencing (WGS) have become more and more affordable. We have seen rapid accumulation of these sequencing data from cohorts with different traits and diseases and are about to see even more in the near future. Now, the major challenge is what we could learn from these data. Specifically, **could we use mutation data to learn what types of mutations are relevant to disease eitiology and which genes are likely to be risk genes**. TADA-A is a statistical framework designed to answer these important questions using *de novo mutations* (DNM), one type of mutations that spontaneously arise in an offspring and not present in its parents. The two main types of input data for TADA-A is DNM and functional/conservational annotation data. It works by integrating DNM information across different studies while accounting for technical variations that might affect observed DNM mutation rates among these studies. It first adjusts baseline mutation rates for each study and estimates the effect sizes of annotations using DNM data from all studies. Finally, it uses relevant annotations to predict the deleteriouss of each mutation and predict disease risk genes. 

## 2. Prerequisites

### 2.1 bedtools
Bedtools need to be installed and added to your PATH. We suggest using v2.17.0, which has been provided in the companion files for you to download. Other versions might be incompatible because of modification of input arguments for some of the sub-functions of bedtools. 

### 2.2 bigWigAverageOverBed
This executable has been added in the `external_tools` folder.

## 3. User guide

### 3.1 Step 1: Adjust mutation rates for each study.
Make a R/Rmd file in the `analysis` folder. And follow the examples below to build your own analysis pipeline.

```r
source("../lib/TADA_Annotation.R")
TADA_A_adjust_mutation_rate(mut_file = "../data/Yuen_NM2015_cases_DNM.bed",
                            window_file = "../data/Example_windows_with_div_score.bed",
                            sample_size = 162, 
                            scale_features = c("GC_content", "div_score"),
                            scaling_file_name = "../data/Example_windows_mutrate_with_div_score_scaling_file_for_Yuen_NM2015_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../data/Kong_cases_DNM.bed",
                            window_file = "../data/Example_windows_with_div_score.bed",
                            sample_size = 78, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../data/Example_windows_mutrate_with_div_score_scaling_file_for_Kong_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../data/Wu_cases_DNM.bed",
                            window_file = "../data/Example_windows_with_div_score.bed",
                            sample_size = 32, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../data/Example_windows_mutrate_with_div_score_scaling_file_for_Wu_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../data/Jiang_cases_DNM.bed",
                            window_file = "../data/Example_windows_with_div_score.bed",
                            sample_size = 32, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../data/Example_windows_mutrate_with_div_score_scaling_file_for_Jiang_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../data/Michaelson_cases_DNM.bed",
                            window_file = "../data/Example_windows_with_div_score.bed",
                            sample_size = 10, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../data/Example_windows_mutrate_with_div_score_scaling_file_for_Michaelson_cases_DNM.txt")
```

In the example code above, we used `TADA_A_adjust_mutation_rate` to adjust the baseline observed mutation rate for each study. We used DNM data from five different studies.

The documentations of the parameters of `TADA_A_adjust_mutation_rate` are listed below

```
--mut_file

A DNM bed file with three columns, separated by "\t"
chr \t start(0-based) \t end(1-based)

--window_file

A window file that has several columns. Each row is a 50-bp window with coordinates information, its associated gene name and other features that are used to adjust baseline mutation rates.


```

|chr	|start	|end	|site_index	|genename	|mutrate	|coding	|promoter	|GC_content	|div_score|
|----|----|----|----|----|----|----|----|----|----|
|chr1	|68090	|68140	|1	|OR4F5	|7.29662e-07	|0	|1	|0.38	|0.0740137843475107|
|chr1	|68140	|68190	|2	|OR4F5	|5.25763e-07	|0	|1	|0.38	|0.0740137843475107|
|chr1	|68190	|68240	|3	|OR4F5	|7.657e-07	|0	|1	|0.42	|0.0740137843475107|


