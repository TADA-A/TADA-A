# a_temp <- new_control$data_matrix[new_control$data_matrix$Prediction == "dnv_proband" & new_control$data_matrix$Type == "SNV" & is.element(new_control$data_matrix$ID, study_ID),1:5]
# a_temp[,2] <- as.numeric(as.character(a_temp[,2]))
# a_temp[,2] <- a_temp[,2] - 1
# write.table(a_temp, "/media/yuwen/F/TADA-A/data/Michaelson_cases_DNM_with_allele_info.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
# 
# 
# 
# 
# mut_files = c("../data/Yuen_NM2015_cases_DNM_with_allele_info.txt","../data/Kong_cases_DNM_with_allele_info.txt","../data/Wu_cases_DNM_with_allele_info.txt", "../data/Jiang_cases_DNM_with_allele_info.txt", "../data/Michaelson_cases_DNM_with_allele_info.txt")
# window_file = "../data/Example_windows.bed"
# mutrate_scaling_files = c("../data/Example_windows_mutrate_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../data/Example_windows_mutrate_scaling_file_for_Kong_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Wu_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Jiang_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Michaelson_cases_DNM.txt")
# sample_sizes = c(162, 78, 32, 32, 10)
# gene_prior_file = "../data/Example_gene_prior.txt"
# nonAS_noncoding_annotations = c("../data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed","../data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed")
# AS_noncoding_annotations = list(c("../data/spidex_public_noncommercial_v1_0.tab_alt_A_lower10pct.bed", "../data/spidex_public_noncommercial_v1_0.tab_alt_C_lower10pct.bed", "../data/spidex_public_noncommercial_v1_0.tab_alt_G_lower10pct.bed","../data/spidex_public_noncommercial_v1_0.tab_alt_T_lower10pct.bed"))
# report_proportion = 100/18665
# chunk_partition_num =1
# node_n = 6
# mutrate_ref_files = c("../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
#                       "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
#                       "../other_annotatio ns/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
#                       "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw")


# define a function that is used in the allele-specific model
# 
TADA_A_read_info <- function(mut_files = c("../data/Yuen_NM2015_cases_DNM_with_allele_info.txt","../data/Kong_cases_DNM_with_allele_info.txt","../data/Wu_cases_DNM_with_allele_info.txt", "../data/Jiang_cases_DNM_with_allele_info.txt", "../data/Michaelson_cases_DNM_with_allele_info.txt"),
                             window_file = "../data/Example_windows.bed",
                             mutrate_scaling_files = c("../data/Example_windows_mutrate_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../data/Example_windows_mutrate_scaling_file_for_Kong_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Wu_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Jiang_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Michaelson_cases_DNM.txt"),
                             sample_sizes = c(162, 78, 32, 32, 10),
                             gene_prior_file = "../data/Example_gene_prior.txt",
                             nonAS_noncoding_annotations = c("../data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed","../data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed"),
                             AS_noncoding_annotations = NA,
                             report_proportion = 100/18665,
                             chunk_partition_num =1 ,
                             node_n = 6,
                             mutrate_ref_files = c("../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw")){

# [mut_file] is a vector of files with DNM infomation in a txt format. The first three columns are chromosome, 0-based start and 1-based end, followed by two columns of ref and alt alleles.
# The code currently only works for SNVs. 
# The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
# [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end" followed by features that might affect local background mutation rate.
# [sample_sizes] is a vector of the number of individuals in each study. The order of the numbers should match the order of mutation files in [mut_files]
# [mutrate_scaling_files] is a vector of files that have the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. 
# each mutrate_scaling_file matches to one pair of window_file and mut_file. The order of files in [mutrate_scaling_file] should match that in [mut_files]
# [nonAS_noncoding_annotations] a vector of non-allele-specific non-coding annotations. Each element is a name of the file that has one non-coding annotation in BED format. Non-coding annotations that are not overlapped with regions in [window_file] will not be used in model fitting. 
# [AS_noncoding_annotations] ia NA or a list of vectors of allele specific annotations. i.e., Each type of noncoding annotation is an element in the list. An element is comprised of 4 different bed files, corresponding to noncoding annotatins of 
# this type based on the alternative allele, A, T, C, or G. e.g., "spidex_lower10pct_alt_A.bed" is a bed file that has genomic intervals representing the union of all bases which, if mutated to an A allele, have a spidex score lower than 10pct of all
# possible spidex scores. 
# [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. 
# [report_proportion] Choose the top X% TADA genes to estimate RR. 
# [node_n] is the number of nodes used to run a certain chunk of the code, default is 6
# [mutrate_ref_files] is a vector of mutrate files in the bigwiggle format. These files have base-level mutation rates to a specific allele, A, C, G, T. 

  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  #mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- windows
  # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  for(i in 1:length(mut_files)){
    mutrate_scaling <- fread(mutrate_scaling_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    system(paste("echo \"Finished reading mutrate scaling file ", mutrate_scaling_files[i], ".\n\"", sep = ""))
    system("date")
    coverage <- coverage[mutrate_scaling, on = "site_index"]
    colnames(coverage)[length(colnames(coverage))] <- paste("scaling_factor_study_", i, sep = "")
  }
  
  # get the piror probability of genes.
  gene_prior <- fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # merge gene prior info
  coverage <- coverage[gene_prior, on = "genename"]
  coverage <-coverage[complete.cases(coverage)]
  
  # select genes based on TADA prior probability and [report_proportion]
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(nrow(genes_for_report)*report_proportion)]
    coverage <- coverage[genes_for_report, on = "genename"]
    coverage <-coverage[complete.cases(coverage)]
  }else{
    genes_for_report  <- gene_prior[,1] # choose all the genes in TADA coding table 
  }
  
  # now need to extropolate window mutation file to base level mutation file
  coverage <-coverage[,c("chr","start","end",paste("scaling_factor_study_", seq(1,length(mut_files)), sep = ""),"genename"),with = FALSE]
  coverage$ID <- paste(coverage$genename, coverage$start, sep = "_")
  
  total_rows <- nrow(coverage)
  interval <- floor(total_rows/chunk_partition_num)
  data_bins <- c(rep(seq(1,chunk_partition_num), each = interval),rep(chunk_partition_num, total_rows -interval*chunk_partition_num))
  
  # split into 20 different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long
  coverage <- split(coverage, data_bins)
  
  #funtion to expand windows to bases
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row["genename"],start,sep = "_"), table_row["ID"])
  }
  
  options(warn=-1)
  cl <- makeCluster(node_n)
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  
  # get nonAS feature number
  nonAS_feature_number <- length(nonAS_noncoding_annotations)
  # get AS feature number
  AS_feature_number <- length(AS_noncoding_annotations)
  
  # get total feature number
  feature_number = nonAS_feature_number + AS_feature_number
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,4:(4 + feature_number - 1)],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,4:(4 + feature_number - 1)])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  
  # build a list to store data
  data_partition <-list()
  alt_letters <- c("A","C","G","T")
  # here only deals with the noncoding parts that are within 10kb of TSSs of genes. Will deal with 
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage[[i]], 1, window_expansion))
    system(paste("echo \"Finished partitioning base-level coordinates data at Round ", i, ".\n\"", sep = ""))
    system("date")
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
    
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    # read in allele-specific base-level mutation rate
    for(j in 1:length(mutrate_ref_files)){
      command <- paste("../external_tools/bigWigAverageOverBed ", mutrate_ref_files[j], " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
      system(command)
      command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
      system(command)
      coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
      colnames(coverage_noncoding_base_mutrate) <- c("base_ID",paste("base_mutrate_alt_", alt_letters[j], sep = ""))
      coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[coverage_noncoding_base_mutrate, on = "base_ID"]
      system(paste("echo \"Finished obtaining base-level mutation rate for alt allele ", alt_letters[j], ".\n\"", sep = ""))
      system("date")
    }
    
    # read in non allele-specific epigenomic annotations
    epi_ID = 1
    if (!is.na(nonAS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in nonAS_noncoding_annotations){
        command <- paste("bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, sep = "_"))
        coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        system(paste("echo \"Finished reading non-allele specific noncoding annotations ", epi_ID, ".\n\"", sep = ""))
        system("date")
        epi_ID <- epi_ID + 1
      }
    }
    
    # read in allele-specific epigenomic annotations, now the base assignment for such noncoding annotations is based on the 50-bp windows. 
    # Should not be a big problem as distal introns are assigned based on refseq_ID for intron regions, should be consistent with assignment by spidex.
    # though there is a small chance that a gene's distal intron, which is close to the promoter of another gene and within 10 kb from the TSSs of the latter gene, may be mis-assigned to the latter gene.
    # To completely correct for this issue, need to allow epigenetic annotation to take its own gene assignment, which might be necessary in some situations under strict criteria, such as splicing annotaion. 
    
    if (!is.na(AS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in AS_noncoding_annotations){
          for(k in 1:length(epi)){
            command <- paste("bedtools coverage -a ", epi[k], " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
            system(command)
            base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
            base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
            colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, alt_letters[k], sep = "_"))
            coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
          }
        system(paste("echo \"Finished reading allele specific noncoding annotations ", epi_ID, ".\n\"", sep = ""))
        system("date")
        epi_ID <- epi_ID + 1
      }
    }
    
   
    # now for each study, read in data, and collapse data based on noncoding annotation configuration
    for(m in 1:length(mut_files)){
      mut <- fread(mut_files[m], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      for(letter in alt_letters){
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        # a very important step here is removing bases that have adjusted_mutrate_base 0. This happens when the allele of the mutrate we are using is just the reference allele. 
        # By doing this, we make the computation of likelihood valid, also we automatically removed bases with nonAS annotations but with mutant allele the same with ref allele at this step
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage_noncoding_for_base_mutrate_temp[[paste("base_mutrate_alt_", letter, sep = "")]] !=0]
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage[[i]][,c(paste("scaling_factor_study_", m, sep = ""), "genename","ID"), with = FALSE],on = "ID"]
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- coverage_noncoding_for_base_mutrate_temp[, paste("base_mutrate_alt_", letter, sep = ""), with = FALSE] * coverage_noncoding_for_base_mutrate_temp[, paste("scaling_factor_study_", m, sep = ""), with = FALSE] * 2 * sample_sizes[m]
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        mut_allele <- mut[V5 == letter]
        fwrite(mut_allele, paste(prefix, "_temp_mut_allele.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        command <- paste("bedtools coverage -a ", paste(prefix, "_temp_mut_allele.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""),sep = "")
        system(command)
        base_with_mut <- fread(paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_with_mut <- base_with_mut[,c("V4","V5"), with = FALSE]
        colnames(base_with_mut) <- c("base_ID", "mut_count")
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[base_with_mut, on = "base_ID"]
        
        # have to collpase data at this point, otherwise, the I will run out of RAM if process 1000 genes every time. 
        
        # first partition by gene for the current chunk
        coverage_noncoding_for_base_mutrate_temp<- split(coverage_noncoding_for_base_mutrate_temp, coverage_noncoding_for_base_mutrate_temp$genename)
        # then partition by feature configuration for each gene in the current chunk
        coverage_noncoding_for_base_mutrate_temp <- sapply(coverage_noncoding_for_base_mutrate_temp, partition_feature, simplify = FALSE)
        # add compact data
        data_partition <- append(data_partition, coverage_noncoding_for_base_mutrate_temp)
        system(paste("echo \"Finished read in mutation data and make them into the compact format for Study ", m, " and allele ", letter, ".\n\"", sep = ""))
        system("date")
      }
    }
  }
  
  stopCluster(cl)
  options(warn = 0)
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste("echo \"Temp files cleaned and data recording finished!\n\""))
  system("date")
  return(list(base_info = data_partition))
}


TADA_A_RR_estimate <-function(data, selected_annotations, gene_prior_file, optimization_iteration = 2000){
  #[data] is the [base_info] returned from [TADA_A_reading_in_annotations], which contains all the allele specific data across all studies
  #[selected_annotations] is a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
  #[gene_prior_file], #[gene_prior_file], a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  gene_prior$prior <- as.numeric(gene_prior$prior)
  gene_prior$genename = as.character(gene_prior$genename)
  # notice the fr function is different from the function that deals with dataset without partition, needs to consider splicing mutation together
  fr <-function(x){ # x is the RRs of selected noncoding annotations. 
    all_rr = x # all_rr is the all regulatory features not including splicing
    cal_logP_Zg1 <- function(data_partition_element){
      cal_logP_Zg1_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[2]]+data_partition_element_level2[[4]]*data_partition_element_level2[[1]][selected_annotations]%*%all_rr-data_partition_element_level2[[3]]*exp(data_partition_element_level2[[1]][selected_annotations]%*%all_rr)-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, cal_logP_Zg1_level2))
    }
    
    cal_logP_Zg0 <- function(data_partition_element){
      cal_logP_Zg0_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[2]]-data_partition_element_level2[[3]]-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, cal_logP_Zg0_level2))
    }
    
    logP_Zg1 = sapply(data, cal_logP_Zg1)
    logP_Zg0 = sapply(data, cal_logP_Zg0)
      
    logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(logP_Zg1))
    logP_table <- logP_table[gene_prior, on = "genename"]
    logP_table <- logP_table[complete.cases(logP_table)]
    ll_sum1 <- sum(by(logP_table, logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
    ll_sum1
  }
  # Use [optim] to do optimization, for non-splicing_mutations
  feature_number <- length(selected_annotations)
  if (feature_number == 1){
    mle <- optim(rep(0.1, feature_number), fr, method = "Brent", lower = -1, upper = 10,control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = TRUE)
  }else{
    mle <- optim(rep(0.1, feature_number), fr ,control=list("fnscale"=-1, "maxit" = optimization_iteration), hessian = TRUE)
  }
  
  # get confidence intervals of RR estimates
  fisher_info <- solve(-mle$hessian)
  prop_sigma <- sqrt(c(diag(fisher_info)))
  rr_estimate <- c(mle$par)
  upper<-rr_estimate+1.96*prop_sigma
  lower<-rr_estimate-1.96*prop_sigma
  rr_report <-data.frame(relative_risk = rr_estimate, lower_bound = lower, upper_bournd = upper)
  
  list(mle = mle, rr_report = rr_report)
}
