# TADA-A
A statistical framework for mapping risk genes from *de novo* mutations in whole-genome sequencing studies

## 1. Introduction
With the fast pace of technology revolution in the field of genomics, whole-exome sequencing (WES) and whole-genome sequencing (WGS) have become more and more affordable. We have seen rapid accumulation of these sequencing data from cohorts with different traits and diseases and are about to see even more in the near future. Now, the major challenge is what we could learn from these data. Specifically, **could we use mutation data to learn what types of mutations are relevant to disease eitiology and which genes are likely to be risk genes**. TADA-A is a statistical framework designed to answer these important questions using *de novo mutations* (DNM), one type of mutations that spontaneously arise in an offspring and not present in its parents. The two main types of input data for TADA-A is DNM and functional/conservational annotation data. It works by integrating DNM information across different studies while accounting for technical variations that might affect observed DNM mutation rates among these studies. It first adjusts baseline mutation rates for each study and estimates the effect sizes of annotations using DNM data from all studies. Finally, it uses relevant annotations to predict the deleteriouss of each mutation and predict disease risk genes. 

## 2. Prerequisites

### 2.1 bedtools
We used bedtools v2.17.0, which has been provided with TADA-A. Other versions might be incompatible because of modification of input arguments for some of the sub-functions of bedtools. But you don't need to worry about this issue. Executables are included in `external_tools/bedtools-2.17.0/bin/` and ready to use after loading the R functions. 

### 2.2 bigWigAverageOverBed
This executable has been added in the `external_tools` folder. 

### 2.3 R parckages
`data.table` and `parallel` need to be installed in R. 

### 2.4 other mutation and annotation files
You need to download and extract `MS_data.tar.gz` in the `TADA-A` folder to replicate the major analyses performed in the manuscript. The link is https://drive.google.com/open?id=13_ofZtsco5SCsqgwdMvZQr_BwBLxffRu. 
```
# after downloading, put MS_data.tar.gz in the main folder `TADA-A`
tar -xvzf MS_data.tar.gz
```
You need to have 10GB space to hold the uncompressed files. 

## 3. Test run

We provided a `test_run.Rmd` in the `test_data` folder. Running the `test_run.Rmd` would generate a `test_run.html`, which should be the same as the `test_run_results.html` in the same folder. It takes about 5mins to finish (Intel(R) Core(TM) i7-4790K CPU @ 4.00GHz).

## 4. User guide
Follow this guide, you will get familiar with the usage of TADA-A and replicate the major analyses done in the accompanying manuscript.

**Note:**
**To analyze data, we suggest put all R or Rmd files with the executing code below in the `analysis` folder under `TADA-A`. Invoke R in the `analysis` folder.**

### 4.1 Step 1: Adjust mutation rates for each study.
Make a R/Rmd file in the `analysis` folder. And follow the examples below to build your own analysis pipeline.

```r
library(data.table)
library(parallel)
source("../lib/TADA_Annotation.R")
TADA_A_adjust_mutation_rate(mut_file = "../MS_data/Yuen_NM2015_cases_DNM.bed",
                            window_file = "../MS_data/Example_windows_with_div_score.bed",
                            sample_size = 162, 
                            scale_features = c("GC_content", "div_score"),
                            scaling_file_name = "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Yuen_NM2015_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../MS_data/Kong_cases_DNM.bed",
                            window_file = "../MS_data/Example_windows_with_div_score.bed",
                            sample_size = 78, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Kong_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../MS_data/Wu_cases_DNM.bed",
                            window_file = "../MS_data/Example_windows_with_div_score.bed",
                            sample_size = 32, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Wu_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../MS_data/Jiang_cases_DNM.bed",
                            window_file = "../MS_data/Example_windows_with_div_score.bed",
                            sample_size = 32, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Jiang_cases_DNM.txt")

TADA_A_adjust_mutation_rate(mut_file = "../MS_data/Michaelson_cases_DNM.bed",
                            window_file = "../MS_data/Example_windows_with_div_score.bed",
                            sample_size = 10, 
                            scale_features = c("GC_content", "div_score"), 
                            scaling_file_name = "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Michaelson_cases_DNM.txt")
```

In the example code above, we used `TADA_A_adjust_mutation_rate` to adjust the baseline observed mutation rate for each study. We used DNM data from five different studies.

The documentations of the parameters of `TADA_A_adjust_mutation_rate` are listed below


--mut_file

A string representing a DNM bed file with three columns, separated by "\t". Notice, here the mutation file does not need to have allele information.
chr \t start(0-based) \t end(1-based)

--window_file

A window file that has several columns. Each row is a 50-bp window with coordinates information, uncalibrated baseline mutation rates, its associated gene name and other features that are used to adjust baseline mutation rates. A example window file is provided as `../MS_data/Example_windows_with_div_score.bed` and the first four rows are shown below. If you use your own window file,  the first 6 columns are required and the names of these columns should not be changed.  The rest of the columns are features that are used to adjust mutation rates so could be customized. 


|chr	|start	|end	|site_index	|genename	|mutrate	|coding	|promoter	|GC_content	|div_score|
  |----|----|----|----|----|----|----|----|----|----|
  |chr1	|68090	|68140	|1	|OR4F5	|7.29662e-07	|0	|1	|0.38	|0.0740137843475107|
  |chr1	|68140	|68190	|2	|OR4F5	|5.25763e-07	|0	|1	|0.38	|0.0740137843475107|
  |chr1	|68190	|68240	|3	|OR4F5	|7.657e-07	|0	|1	|0.42	|0.0740137843475107|
  
  --sample_size

A integer representing the number of individuals in the study

--scale_features

A vector showing the names of the features for which normalization are needed when adjusting mutation rates. 

--scaling_file_name

A string giving the name of the output file that gives a mutation rate scaling factor for each window, which will be used in the following steps. The first three rows of an example output file is shown below. As you can see, we use `site_index` as the main identifier to link the window file and the mutation rate scaling file. 

|site_index	|scaling_factor|
  |----|----|
  |1	|0.454254325330186|
  |2	|0.454254325330187|
  
  
  ### 4.2 Step 2: Feature selection of annotations.
  
  #### 4.2.1 Read in DNM data and annotation
  
  We read DNM data and annotation data using `TADA_A_read_info` and store all the information into `compact_data`, a compact data form benefiting from our categorization trick. This categorization trick greatly reduced the size of the data and facilitates fast parameter inference. 
```r
compact_data <- TADA_A_read_info(mut_files = c("../MS_data/Yuen_NM2015_cases_DNM_with_allele_info.txt","../MS_data/Kong_cases_DNM_with_allele_info.txt","../MS_data/Wu_cases_DNM_with_allele_info.txt", "../MS_data/Jiang_cases_DNM_with_allele_info.txt", "../MS_data/Michaelson_cases_DNM_with_allele_info.txt"),
                                 window_file = "../MS_data/Example_windows_with_div_score.bed",
                                 mutrate_scaling_files = c("../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Kong_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Wu_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Jiang_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Michaelson_cases_DNM.txt"),
                                 sample_sizes = c(162, 78, 32, 32, 10),
                                 gene_prior_file = "../MS_data/Example_gene_prior.txt",
                                 nonAS_noncoding_annotations = c("../MS_data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed", "../MS_data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed", "../MS_data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed","../MS_data/170820_gerp_gt2.all.sorted.merged.bed","../MS_data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr_gerp_gt2.bed", "../MS_data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr_gerp_gt2.bed", "../MS_data/Encode_DHS_union_within_10kb_and_promoter_no_utr_gerp_gt2.bed"),
                                 AS_noncoding_annotations = list(c("../MS_data/spidex_public_noncommercial_v1_0.tab_alt_A_lower10pct.bed", "../MS_data/spidex_public_noncommercial_v1_0.tab_alt_C_lower10pct.bed", "../MS_data/spidex_public_noncommercial_v1_0.tab_alt_G_lower10pct.bed","../MS_data/spidex_public_noncommercial_v1_0.tab_alt_T_lower10pct.bed"), c("../MS_data/whole_genome_SNVs_gt15_altA_within_10kb_and_promoter_no_utr.bed", "../MS_data/whole_genome_SNVs_gt15_altC_within_10kb_and_promoter_no_utr.bed", "../MS_data/whole_genome_SNVs_gt15_altG_within_10kb_and_promoter_no_utr.bed", "../MS_data/whole_genome_SNVs_gt15_altT_within_10kb_and_promoter_no_utr.bed"), c("../MS_data/whole_genome_SNVs_gt15_altA_within_10kb_and_promoter_no_utr_with_Noonan_brain_roadmap_union.bed", "../MS_data/whole_genome_SNVs_gt15_altC_within_10kb_and_promoter_no_utr_with_Noonan_brain_roadmap_union.bed","../MS_data/whole_genome_SNVs_gt15_altG_within_10kb_and_promoter_no_utr_with_Noonan_brain_roadmap_union.bed","../MS_data/whole_genome_SNVs_gt15_altT_within_10kb_and_promoter_no_utr_with_Noonan_brain_roadmap_union.bed"), c("../MS_data/whole_genome_SNVs_gt15_altA_within_10kb_and_promoter_no_utr_with_E081_E082_intersection.bed", "../MS_data/whole_genome_SNVs_gt15_altC_within_10kb_and_promoter_no_utr_with_E081_E082_intersection.bed", "../MS_data/whole_genome_SNVs_gt15_altG_within_10kb_and_promoter_no_utr_with_E081_E082_intersection.bed", "../MS_data/whole_genome_SNVs_gt15_altT_within_10kb_and_promoter_no_utr_with_E081_E082_intersection.bed"), c("../MS_data/whole_genome_SNVs_gt15_altA_within_10kb_and_promoter_no_utr_with_Encode_DHS_union.bed", "../MS_data/whole_genome_SNVs_gt15_altC_within_10kb_and_promoter_no_utr_with_Encode_DHS_union.bed", "../MS_data/whole_genome_SNVs_gt15_altG_within_10kb_and_promoter_no_utr_with_Encode_DHS_union.bed", "../MS_data/whole_genome_SNVs_gt15_altT_within_10kb_and_promoter_no_utr_with_Encode_DHS_union.bed")),
                                 report_proportion = 1000/18665,
                                 chunk_partition_num =1,
                                 node_n = 6,
                                 mutrate_ref_files = c("../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                                                       "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                                                       "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                                                       "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw")
)
```
The documentations of the parameters of `TADA_A_read_info` are listed below

--mut_files

A vector with names of DNM files. Notice here allelic information needs to be included as our model is allele-aware. Below is an example of the first three rows of one DNM file.

|    |    |    |    |    |
  |----|----|----|----|----|
  |chr3	|26636920	|26636921	|A	|G|
  |chr3	|70569020	|70569021	|A	|C|
  |chr3	|82677644	|82677645	|A	|G|
  
  --window_file

Must be the same window file used in the mutation rate adjustion step (Step 3.1).

--mutrate_scaling_files

A vector with the names of mutation rate scaling files which were generated in the mutation rate adjustion step (Step 3.1). The ordering of these files should be consistent with that of the mutation rate files specified by `--mut_files`

--sample_sizes

A vector with the sample sizes of all the studies. The ordering of these numbers should be consistent with the mutation rate files specified by `--mut_files`

--gene_prior_file

A string representing the name of a file with the prior probability of each gene being risk gene. Below is an example of the first three rows of one prior file. 

|    |    | 
  |----|----|
  |genename	|prior|
  |CHD8	|0.999999999694329|
  |SCN2A	|0.999999999591982|
  
  --nonAS_noncoding_annotations

A vector representing the names of non-allele specific annotations. All annotations should be in a 3-column BED format. Each annotation file covers all the bases (or at least all the bases included in the window file) that have that annotation.

--AS_noncoding_annotations

A list representing allele-specific annotations. Each element of the list is a four-element vector of the names of allele-specific annotation files for each annotation. Each allele-specific annotation file of an annotation is in 3-column BED format, covering all the bases (or at least all the bases included in the window file) that have that annotation if mutated to a certain allele. The ordering of these allele-specific annotation files for one annotation strictly follows mutant allele as "A", "C", "G", and "T".

--report_proportion

A number showing the proportion of genes whose DNM and annotation information will be collected and used in relative risk estimation of annotations. All the genes in the `gene_prior_file` will be ranked from high to low based on their prior probabilities. Only the top proportion based on `report_proportion` of all genes will be used. 

--chunk_partition_num

A number specifying how many partitions will the `window_file` be split into. We do partitions in order to reduce the computation burden by only handling a relatively small set of genomic regions when counting DNM and recording annotations at each base. If only look at top 1000 genes, setting `--chunk_partition_num` to 1 would be usually sufficient.

--node_n

The number of computational nodes that needs to be used.

--mutrate_ref_files

A vector specifying the names of the base-level allele-specific baseline mutation rate files. These files are in the bigWiggle format, storing base-level mutation rates for mutant allele as "A", "C", "G", and "T", sequentially. Users could use their own mutation rate models but need to build their own base-level allele-specific mutation rate files accordingly. 

#### 4.2.2 Feature selection.

Estimate the relative risks of each individual annotation supplied to `TADA_A_read_info`, and only use the sinificant ones for the next round of joint estimation. The code below performs relative risk estimation. We have 13 features here, so we estimate the relative risk for each feature sequentially.

```r
for(i in 1:12){
  print(paste("Relative risk estimate of the ", i, "th annotation"), sep = "")
  print(TADA_A_RR_estimate(data = compact_data$base_info, selected_annotations = c(i), gene_prior_file = "../MS_data/Example_gene_prior.txt", optimization_iteration = 2000))
}
```
The documentations of the parameters of `TADA_A_RR_estimate` are listed below

--data

The `base_info` object returned by `TADA_A_read_info`, which stores DNM information, mutation rates, and annotations for each gene in a compact format. 

--selected_annotations

A vector indicating which annotations are used in relative risk estimation. In the feature selection step, we always estimate relative risk for each feature separately, so the vector only has one number specifying which annotation will be used. As you may member, we have multiple annotations specified by `--nonAS_noncoding_annotations` and `AS_noncoding_annotations`. We numbered these annotations from 1 to the total number of annotations. For example, if we have three nonallele-specific annotations `c("nonAS_Annotation_1", "nonAS_Annotation_2", "nonAS_Annotation_3")` specified by `--nonAS_noncoding_annotations`, and two allele-specific annotations `list(c("AS_Annotation_1_alt_A", "AS_Annotation_1_alt_A", "AS_Annotation_1_alt_A", "AS_Annotation_1_alt_A"), c("AS_Annotation_2_alt_A", "AS_Annotation_2_alt_A", "AS_Annotation_2_alt_A", "AS_Annotation_2_alt_A"))` specified by `AS_noncoding_annotations`. Then we have five annotations together, "nonAS_Annotation_1", "nonAS_Annotation_2", "nonAS_Annotation_3", "AS_Annotation_1", and "AS_Annotation_2" will be numbered 1, 2, 3, 4, 5, respectively. If we want to estimate the relative risk of "nonAS_Annotation_3", we set `--selected_annotations` to `c(3)`. If we want to jointly estimate the relative risks of "nonAS_Annotation_3" and "AS_Annotation_1", we set `--selected_annotations` to `c(3,4)`.

--gene_prior_file

A string representing the name of a file with the prior probability of each gene being risk gene. Below is an example of the first three rows of one prior file. In most scenarios, you want to be consistent with `TADA_A_read_info` regarding to the choice of priors.

|    |    | 
  |----|----|
  |genename	|prior|
  |CHD8	|0.999999999694329|
  |SCN2A	|0.999999999591982|
  
  --optimization_iteration

The maximum number of iterations when performing optimization. `Optim()` was used for optimization. The search space for when only one parameter is estimated is from -1 to 10.

--mode 

A string that is `"regular"` (default), or `"single_fast"`. `"single_fast"` is used when estimating RR from only one annotation ( when running `TADA_A_read_info`, only one annotation is provided) of lots of genes (e.g., all genes), would be at least 5 times faster.

[output]:
  
  The output of `TADA_A_RR_estimate` has two objects. The first one `mle` is the output from `optim()`, the second one `rr_report` is a data.frame. Each row of the data.frame is an annotation. Columns are log(Relative risk), 95% confidence interval lower bound of log(RR), and 95% confidence interval upper bound of log(RR).

#### 4.2.3 Jointly estimating the relative risks of annotations.

We perform joint estimation for annotations that pass the feature selection step. Example code is shown below. We set `--selected_annotations` to be `c(1,5,8)` as the 1st, 5th and 8th annotations passed the feature selection. 

```r
TADA_A_RR_estimate(data = compact_data$base_info, selected_annotations = c(1,5,8), gene_prior_file = "../MS_data/Example_gene_prior.txt", optimization_iteration = 2000)
```


#### 4.2.4 Predict risk genes
Use the relative risks of annotations from Step 3.2.3 to identify risk genes. The example code is below. Notice, if previously in relative risk estimation, we only used top genes (e.g., top 1000 genes based on priors), then `compact_data$base_info`, would only contain information for these top 1000 genes. To predict risk genes for all potential genes, we need to first run `TADA_A_read_info` again over all the genes in the `window_file`. 

The difference compared to Step 3.2.2 is that 1) now we only take in annotations that have passed the feature selection when running `TADA_A_read_info`; 2) set `--report_proportion` to 1 so all genes in the `window_file` would be included; 3) set `--chunk_partition_num` to a bigger number such as `20` to avoid memory overflow issue, which would arise when each chunk has to cover too many genomic positions. 

```r
compact_data_after_feature_selection <- TADA_A_read_info(mut_files = c("../MS_data/Yuen_NM2015_cases_DNM_with_allele_info.txt","../MS_data/Kong_cases_DNM_with_allele_info.txt","../MS_data/Wu_cases_DNM_with_allele_info.txt", "../MS_data/Jiang_cases_DNM_with_allele_info.txt", "../MS_data/Michaelson_cases_DNM_with_allele_info.txt"),
                                                         window_file = "../MS_data/Example_windows_with_div_score.bed",
                                                         mutrate_scaling_files = c("../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Kong_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Wu_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Jiang_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Michaelson_cases_DNM.txt"),
                                                         sample_sizes = c(162, 78, 32, 32, 10),
                                                         gene_prior_file = "../MS_data/Example_gene_prior.txt",
                                                         nonAS_noncoding_annotations = c("../MS_data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../other_annotations/conservation/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr_gerp_gt2.bed"),
                                                         AS_noncoding_annotations = list(c("../MS_data/spidex_public_noncommercial_v1_0.tab_alt_A_lower10pct.bed", "../MS_data/spidex_public_noncommercial_v1_0.tab_alt_C_lower10pct.bed", "../MS_data/spidex_public_noncommercial_v1_0.tab_alt_G_lower10pct.bed","../MS_data/spidex_public_noncommercial_v1_0.tab_alt_T_lower10pct.bed")),
                                                         report_proportion = 18665/18665,
                                                         chunk_partition_num =20,
                                                         node_n = 6,
                                                         mutrate_ref_files = c("../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                                                                               "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                                                                               "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                                                                               "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw")
)
```

With the information of all the genes recorded, we run `TADA_A_get_BFs` to get the risk gene prediction table. Notice now, the numbering of the three anntations become 1, 2, and 3, respectively. This is because we only used these three mutations when running `TADA_A_read_info`. We specify the relative risks in the log scale to `rr` for these three annotations sequentially in a vector. We also provide a Bayes factor table based on independent external information to help boost the power of risk gene prediction. 

```r
TADA_A_get_BFs(data = compact_data_after_feature_selection$base_info, selected_annotations = c(1,2,3), rr = c(0.43,0.80, 1.17), additional_BF_file = "../MS_data/Example_gene_coding_BF.txt")
```

**Notice**
  
  Running `TADA_A_read_info` for all the genes over many annotations may be memory intensive and takes a couple of days. We do provide an alternative strategy here. We could first partition the window file specified by `--window_file` in `TADA_A_read_info` into 50 parts (We provided the partitions in `MS_data` folder as `Example_windows_with_div_score_no_header.bed.temp.partition*.with_header.txt`). For each partition, we run simultaneously on a server `TADA_A_read_info` and store the output `base_info` in `RDS` format. It takes about 2hrs if using ~10 annotations for DNMs from five different studies (with less than 200 individuals each). We then read in the `.RDS` file iteratively and `append()` them together, and use the final object as the input for `--data` for `TADA_A_get_BFs` . Below is an example code.
  

```r
# Read in data from one genomic partition, and save as a .RDS file
# This is just an example code to show how to run analysis for one genomic partition.
# In practice, we could analyze all the genomic partitions the same way simultaneously on a cluster. 

compact_data_00 <- TADA_A_read_info(mut_files = c("../MS_data/Yuen_NM2015_cases_DNM_with_allele_info.txt","../MS_data/Kong_cases_DNM_with_allele_info.txt","../MS_data/Wu_cases_DNM_with_allele_info.txt", "../MS_data/Jiang_cases_DNM_with_allele_info.txt", "../MS_data/Michaelson_cases_DNM_with_allele_info.txt"),
                                   window_file = "../MS_data/windows_partition/Example_windows_with_div_score_no_header.bed.temp.partition00.with_header.txt",
                                   mutrate_scaling_files = c("../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Kong_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Wu_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Jiang_cases_DNM.txt", "../MS_data/Example_windows_mutrate_with_div_score_scaling_file_for_Michaelson_cases_DNM.txt"),
                                   sample_sizes = c(162, 78, 32, 32, 10),
                                   gene_prior_file = "../MS_data/Example_gene_prior.txt",
                                   nonAS_noncoding_annotations = c("../MS_data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed",  "../other_annotations/conservation/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr_gerp_gt2.bed"),
                                   AS_noncoding_annotations = list(c("../MS_data/spidex_public_noncommercial_v1_0.tab_alt_A_lower10pct.bed", "../MS_data/spidex_public_noncommercial_v1_0.tab_alt_C_lower10pct.bed", "../MS_data/spidex_public_noncommercial_v1_0.tab_alt_G_lower10pct.bed","../MS_data/spidex_public_noncommercial_v1_0.tab_alt_T_lower10pct.bed")),
                                   report_proportion = 18665/18665,
                                   chunk_partition_num =1,
                                   node_n = 1,
                                   mutrate_ref_files = c("../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                                                         "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                                                         "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                                                         "../MS_data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw"),
                                   MPI = 00
)

saveRDS(compact_data_00,"../MS_data/genomic_partition_00.RDS")

```



```r
# append .RDS files from all partitions together
data_partition <- list()
test = sprintf("%02d", 0:49)
for(i in 1:50){
  temp <- readRDS(paste("../MS_data/genomic_partition_", test[i], ".RDS", sep = ""))
  data_partition <- append(data_partition, temp$base_info)
}

table <- TADA_A_get_BFs(data = data_partition, selected_annotations = c(1,2,3), rr = c(0.43,0.80, 1.17), additional_BF_file = "../MS_data/Example_gene_coding_BF.txt")
```


The documentations of the parameters of `TADA_A_get_BFs` are listed below

--data

The `base_info` object returned from `TADA_A_RR_estimate`, which stores DNM information, mutation rates, and annotations for each gene in a compact format. 


--selected_annotations

A vector indicating which annotations are selected in risk gene prediction. When using `TADA_A_read_info` to record information for all the genes, we have multiple annotations specified by `--nonAS_noncoding_annotations` and `AS_noncoding_annotations`. We numbered these annotations from 1 to the total number of annotations. For example, if we have three nonallele-specific annotations `c("nonAS_Annotation_1", "nonAS_Annotation_2", "nonAS_Annotation_3")` specified by `--nonAS_noncoding_annotations`, and two allele-specific annotations `list(c("AS_Annotation_1_alt_A", "AS_Annotation_1_alt_A", "AS_Annotation_1_alt_A", "AS_Annotation_1_alt_A"), c("AS_Annotation_2_alt_A", "AS_Annotation_2_alt_A", "AS_Annotation_2_alt_A", "AS_Annotation_2_alt_A"))` specified by `AS_noncoding_annotations`. Then we have five annotations together, "nonAS_Annotation_1", "nonAS_Annotation_2", "nonAS_Annotation_3", "AS_Annotation_1", and "AS_Annotation_2" will be numbered 1, 2, 3, 4, 5, respectively. If we want to predict risk genes using "nonAS_Annotation_3", we set `--selected_annotations` to `c(3)`. If we want to predict risk genes using "nonAS_Annotation_3" and "AS_Annotation_1", we set `--selected_annotations` to `c(3,4)`.

--rr

A vector giving the log(relative risk) of features selected in `--selected_annotations`. The ordering must be consistent between `--rr` and `--selected_annotations`

--additional_BF_file

The name of a file with Bayes factors from independent external information. For example, we may use Bayes factors based on coding mutations from WES studies when our main analyses was performed using noncoding mutations from WGS studies. If we don't have additional bayes factors to use, we could use a file with all the Bayes factors equal to 1 here.  The first three rows of an `--additional_BF_file` is shown below

|genename	|coding_BF|
|----|----|
|CHD8	|7317632288.7186|
|SCN2A	|17086885921.6339|


[Output]:

`$gene_BF_table`:

|genename	|logBF_noncoding	|logBF_coding	|logBF_all	|BF_noncoding	|BF_coding	|BF_all	|FDR_coding	|FDR_all|
|----|----|----|----|----|----|----|----|----|
|APBB1	|2.789694	|2.325413	|5.115108	|16.276045	|10.23091	|166.51873	|0.4502329	|0.0261426|
|NRXN1	|1.079917	|2.936408	|4.016325	|2.944436	|18.84801	|55.49677	|0.1982470	|0.0814069|
|PNPLA7	|1.114084	|2.863029	|3.977114	|3.046777	|17.51450	|53.36279	|0.2273822	|0.0866883|
|TANC2	|1.140735	|2.798141	|3.938877	|3.129068	|16.41411	|51.36087	|0.2697554	|0.0942703|

+ genename: The name of a gene.
+ logBF_noncoding: The log(Bayes Factor) of the gene using WGS noncoding data.
+ logBF_coding: The log(Bayes factor) of the gene using independent external Bayes factor, here from WES coding data.
+ logBF_all: The summation of logBF_noncoding and logBF_coding.
+ BF_noncoding: The Bayes Factor of the gene using WGS noncoding data.
+ BF_coding: The Bayes factor of the gene using independent external Bayes factor, here from WES coding data.
+ BF_all: The multipication of BF_noncoding and BF_coding.
+ FDR_coding: The FDR based on independent external information.
+ FDR_all: The FDR based on combining the independent external information and the data that we analyzed. 
