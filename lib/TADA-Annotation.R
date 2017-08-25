# function to calibrate background mutation rate for each DNM study.
TADA_A_adjust_mutation_rate <- function(mut_file, window_file, sample_size, scale_features, scaling_file_name){
  # [mut_file] is the file with DNM infomation in a BED format. 0-based start and 1-based end. e.g., "../data/Example_windows.bed"
  # The first 4 columns must be chr, start, end, and site_index of genomic windows.
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end", "genename", and "mutrate", followed by features that might affect local background mutation rate.
  # [sample_size] is the number of individuals. 
  # [scale_features] is a vector of the names of features that need to be scaled, this is recommended to apply to continuous features, such as GC content, to make easier the interpretation of the effect sizes of features. 
  # [scaling_file_name] is the name of the file that has the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. 
  
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE)
  prefix <- paste("tmp/", prefix, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  command <- paste("sed -n '1!p' ", window_file, " | awk {'print $1\"\t\"$2\"\t\"$3\"\t\"$4'} > ", paste(prefix, "_temp_windows.bed", sep = ""), sep = "")
  system(command)
  # get the number of SNVs for each window
  command <- paste("bedtools coverage -a ", mut_file, " -b ", paste(prefix, "_temp_windows.bed", sep = ""), " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
  system(command)
  # read in the file with the number of SNVs for each window
  coverage <- fread(paste(prefix,"_temp_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- coverage[,1:5]
  colnames(coverage) <- c("chr","start","end","site_index","mut_count")
  # read in the window file that has feature annotations that might affect mutation rates
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  feature_number <- dim(windows)[2] - 6
  # merge [window] and [coverage] to link mutation count with feature annotations
  coverage <- coverage[windows[,-1:-3], on="site_index"]
 
  coverage <- coverage[mutrate !=0]
  for(i in 1:length(scale_features)){
    coverage[[scale_features[i]]] <- as.vector(scale(coverage[[scale_features[i]]]))
  }
  # write the formula
  f <- paste("mut_count ~ ", paste(colnames(coverage)[8 : (8 + feature_number -1)], collapse = " + "), " + offset(log(2*mutrate*sample_size))", sep = "")
  #fit mutation rate model using a fixed set of features not including histone modification marks. 
  out.offset <- glm(as.formula(f), family = poisson, data = coverage)
  scaling <- data.table(site_index = coverage$site_index, scaling_factor = out.offset$fitted.values / (2 * coverage$mutrate * sample_size))
  fwrite(scaling, scaling_file_name, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
  # remove intermediate files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  # the return value is the effect sizes of feature annotations on background mutation rates. 
  return(summary(out.offset)$coeff)
}


TADA_A_read_info <- function(mut_file = "../data/Yuen_NM2015_cases_DNM.bed",
                                          window_file = "../data/Example_windows.bed",
                                          mutrate_scaling_file = "../data/Example_windows_mutrate_scaling_file_for_Yuen_NM2015_cases_DNM.txt",
                                          sample_size = 162,
                                          gene_prior_file = "../data/Example_gene_prior.txt",
                                          noncoding_annotations = c("../data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed","../data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed"),
                                          report_proportion = 100/18665,
                                          chunk_partition_num =1,
                                          node_n = 6,
                                          mutrate_ref_file = "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"){
  # [mut_file] is the file with DNM infomation in a BED format. 0-based start and 1-based end. e.g., "../data/Example_windows.bed"
  # The first 4 columns must be chr, start, end, and site_index of genomic windows. The rest of the columns are features that might affect baseline background mutation rates and that need to be adjusted for.
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end" followed by features that might affect local background mutation rate.
  # [sample_size] is the number of individuals.
  # [mutrate_scaling_file] is the name of the file that has the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. 
  # each mutrate_scaling_file matches to one pair of window_file and mut_file.
  # [noncoding_annotations] a vector of non-coding annotations. Each element is a name of the file that has one non-coding annotation in BED format. Non-coding annotations that are not overlapped with regions in [window_file] will not be used in model fitting. 
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. 
  # [report_proportion] Choose the top X% TADA genes to estimate RR. 
  # [node_n] is the number of nodes used to run a certain chunk of the code, default is 6
  # [mutrate_ref_file] the file with base level mutation rate reference. e.g., "/media/yuwen/Elements/mutation_rate_raw/MarkDaly/Daly_mutrate.bw"

  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  mutrate_scaling <- fread(mutrate_scaling_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- windows[mutrate_scaling, on = "site_index"]
  
  # get the adjusted mutation rate per base per individual
  coverage$adjusted_mutrate = coverage$mutrate * coverage$scaling_factor
  
  # get the piror probability of genes.
  gene_prior <- fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # merge gene prior info
  coverage_1 <- coverage[gene_prior, on = "genename"]
  coverage_1 <-coverage_1[complete.cases(coverage_1)]
  
  # select genes based on TADA prior probability and [report_proportion]
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(nrow(genes_for_report)*report_proportion)]
    coverage_1 <- coverage_1[genes_for_report, on = "genename"]
    coverage_1 <-coverage_1[complete.cases(coverage_1)]
  }else{
    genes_for_report  <- gene_prior[,1] # choose all the genes in TADA coding table 
  }
  
  # now need to extropolate window mutation file to base level mutation file
  coverage_1 <-coverage_1[,c("chr","start","end","scaling_factor","genename"),with = FALSE]
  coverage_1$ID <- paste(coverage_1$genename, coverage_1$start, sep = "_")
  
  total_rows <- nrow(coverage_1)
  interval <- floor(total_rows/chunk_partition_num)
  data_bins <- c(rep(seq(1,chunk_partition_num), each = interval),rep(chunk_partition_num, total_rows -interval*chunk_partition_num))
  
  # split into 20 different chunks, then for each chunk split by genes, then by feature configuration. If split by gene at the first level, implementation would take too long
  coverage_1 <- split(coverage_1, data_bins)
  # get feature number
  feature_number = length(noncoding_annotations)
  # function to get effective information of each element of partition_by_gene
  # These information are those necessary to compute log-likelihood in the optimization function
  partition_feature <- function(pbg){
    # input is one element of the list of partition_by_gene
    pbg_split <- split(pbg, pbg[,5:(5 + feature_number - 1)],drop = TRUE)
    feature_combination_number <- length(pbg_split)
    # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
    info_for_each_feature <- function(feature_set){
      list(feature_vector = c(as.numeric(feature_set[1,5:(5 + feature_number - 1)])), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
    }
    sapply(pbg_split, info_for_each_feature,simplify = FALSE)
  }
  
  
  #funtion to expand windows to base level
  window_expansion <- function(table_row){
    start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
    data.frame(table_row[1], start, start+1, paste(table_row[5],start,sep = "_"), table_row[6])
  }
  
  options(warn=-1)
  cl <- makeCluster(node_n)
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  
  # build vectors and list to store results

  data_partition <-list()
  
  # here only deals with the noncoding parts that are within 10kb of TSSs of genes. Will deal with 
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage_1[[i]], 1, window_expansion))
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate$start <- as.integer(coverage_noncoding_for_base_mutrate$start)
    coverage_noncoding_for_base_mutrate$end <- as.integer(coverage_noncoding_for_base_mutrate$end)
    
    # write out a bed file to get base-level mutation rates
    fwrite(coverage_noncoding_for_base_mutrate[,1:4],paste(prefix, "_temp_for_mutrate.bed", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    command <- paste("../external_tools/bigWigAverageOverBed ", mutrate_ref_file, " ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = "" ), sep = "")
    system(command)
    command <- paste("awk {'print $1\"\t\"$4'} ", paste(prefix, "_temp_for_mutrate.bed.mutrate", sep = ""), " > ", paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), sep = "")
    system(command)
    coverage_noncoding_base_mutrate <-fread(paste(prefix, "_temp_for_mutrate.bed.mutrate.txt", sep = ""), header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(coverage_noncoding_base_mutrate) <- c("base_ID","base_mutrate")
    a_temp <- coverage_noncoding_for_base_mutrate[coverage_noncoding_base_mutrate, on = "base_ID"]
    a_temp <- a_temp[coverage_1[[i]], on = "ID"]
    coverage_noncoding_mutrate_adjusted <- a_temp[,c("base_ID","chr","start","end","genename")]
    coverage_noncoding_mutrate_adjusted$adjusted_base_mutrate <- a_temp[,c("base_mutrate")]*a_temp[,c("scaling_factor")]*2*sample_size # remember need to add back sample_size
    rm(a_temp)
    
    # get the mutation count for each base
    command = paste("bedtools coverage -a ", mut_file, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""), " > ", paste(prefix,"_temp_base_level_coverage.bed", sep = ""), sep = "")
    system(command)
    mutation_overlap_base <-fread(paste(prefix,"_temp_base_level_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutation_overlap_base <- mutation_overlap_base[,c("V4","V5")]
    colnames(mutation_overlap_base) <- c("base_ID","mut_count")
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[mutation_overlap_base, on = "base_ID"]
    rm(mutation_overlap_base)
    
    coverage_noncoding_mutrate_adjusted<- coverage_noncoding_mutrate_adjusted[,c("genename","base_ID","mut_count","adjusted_base_mutrate")]
    
    # overlap with multiple epi feature.
    if (noncoding_annotations != "no"){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      epi_ID = 1
      for(epi in noncoding_annotations){
        command <- paste("bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, sep = "_"))
        coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[base_in_epi, on = "base_ID"]
        epi_ID <- epi_ID + 1
      }
    }
    
    # remove bases without mutation rates
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[adjusted_base_mutrate!=0]
    
    # remove rows with NA
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[complete.cases(coverage_noncoding_mutrate_adjusted)]
    
    anno_count <- rep(0, nrow(coverage_noncoding_mutrate_adjusted))
    for(i in 1:feature_number){
      anno_count <- anno_count + coverage_noncoding_mutrate_adjusted[[(5 + i -1)]]
    }
    
    # remove rows(bases) that don't have any non-coding features
    coverage_noncoding_mutrate_adjusted <- coverage_noncoding_mutrate_adjusted[anno_count >0,]
    #then partition by gene for the current chunk
    coverage_noncoding_mutrate_adjusted<- split(coverage_noncoding_mutrate_adjusted, coverage_noncoding_mutrate_adjusted$genename)
    # then partition by feature configuration for each gene in the current chunk
    coverage_noncoding_mutrate_adjusted <- sapply(coverage_noncoding_mutrate_adjusted, partition_feature, simplify = FALSE)
    data_partition <- append(data_partition, coverage_noncoding_mutrate_adjusted)
  }
  
  stopCluster(cl)
  options(warn = 0)
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  return(list(base_info = data_partition))
}


TADA_A_RR_estimate <-function(data_list, selected_annotations, gene_prior_file, optimization_iteration = 2000){
  #[data_list] is a list of objects [base_info] returned from [TADA_A_reading_in_annotations], each element of the list has compact base-level information for each gene.
  # If you only want to estiamte RR from one single study, just assign data_list = list(one_study$base_info), where one_study is returned from [TADA_A_reading_in_annotations].
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
    ll_sum1 <- 0
    for(data in data_list){
      logP_Zg1 = sapply(data, cal_logP_Zg1)
      logP_Zg0 = sapply(data, cal_logP_Zg0)
      
      logP_table<-data.table(logP_Zg1 = logP_Zg1, logP_Zg0 = logP_Zg0, genename = names(logP_Zg1))
      logP_table <- logP_table[gene_prior, on = "genename"]
      logP_table <- logP_table[complete.cases(logP_table)]
      # sum of ll from regulaotry features not including splicing mutations
      ll_sum1 <- ll_sum1 + sum(by(logP_table, logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
    }
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



# Bayesian FDR control (PMID:19822692, Section2.3)

Bayesian.FDR <- function(BF, pi0, alpha=0.05) {
  # convert BFs to PPA (posterior probability of alternative model)
  # [BF]: a sorted vector of BFs (in decreasing order)
  # [pi0]: the prior probability that the null model is true
  # [alpha]: the FDR target
  # [Return]: the q-value of each BF, and the number of findings with q below alpha. 
  pi <- 1-pi0
  q <- pi*BF/(1-pi+pi*BF) # PPA
  q0 <- 1 - q # posterior probability of null model
  
  # the FDR at each PPA cutoff
  n <- length(BF)
  FDR <- numeric(n)
  for (i in 1:n) FDR[i] <- sum(q0[1:i]) / i 
  
  # the cutoff
  t <- 1
  while (t <= length(q0) & mean(q0[1:t]) <= alpha) { t <- t+1 }
  return (list(FDR=FDR, ND=t))
}


# calculate the Bayes factor for each gene based on informative noncodng annotations. 
TADA_A_get_BFs <- function(data_list,selected_annotations, rr, additional_BF_file, TADA_p0 = 0.94){
  # [data_list] is a list of objects [base_info] returned from [TADA_A_reading_in_annotations], each element of the list has compact base-level information for each gene. 
  # [selected_annotations] determines which features in the [data_list] are going to be used for calculating bayes factors. e.g., if [selected_annotations] is c(1,3). Then the first and third element of feature_vector in [data_list] will be selected.
  # [rr] is the estimated relative risks of [selected_annotations]. RRs are estimated from [TADA_A_RR_estimate_from_multiple_studies]. 
  # RRs are in the log scale
  # [additional_BF_file], a file that has additional Bayes factors that could be uesd along with noncoding-derived BFs to increase the power of risk gene mapping. e.g., "../data/Example_gene_coding_BF.txt"
  # [TADA_p0] the proportion of genes being non-risk genes. Default = 0.94 for in the context of ASD risk gene mapping
  
  additional_BF <- fread(additional_BF_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  additional_BF[,2] = log(additional_BF[,2])
  colnames(additional_BF) = c("genename", "logBF_coding")
  
  cal_logP_Zg1 <- function(data_partition_element){
    cal_logP_Zg1_level2 <-function(data_partition_element_level2){
      data_partition_element_level2[[2]]+data_partition_element_level2[[4]]*data_partition_element_level2[[1]][selected_annotations]%*%rr-data_partition_element_level2[[3]]*exp(data_partition_element_level2[[1]][selected_annotations]%*%rr)-data_partition_element_level2[[5]]
    }
    sum(sapply(data_partition_element, cal_logP_Zg1_level2))
  }
  
  cal_logP_Zg0 <- function(data_partition_element){
    cal_logP_Zg0_level2 <-function(data_partition_element_level2){
      data_partition_element_level2[[2]]-data_partition_element_level2[[3]]-data_partition_element_level2[[5]]
    }
    sum(sapply(data_partition_element, cal_logP_Zg0_level2))
  }
  
  # define a list to hold log_noncoding Bayes factor from each study
  logBF_noncoding_list <- list()
  for(data in data_list){
    logP_Zg1 <- sapply(data, cal_logP_Zg1)
    logP_Zg0 <- sapply(data, cal_logP_Zg0)
    
    # define a function to flag genes that don't have any bases with any informative noncoding annotations that are considered to affect relative risk (i.e., noncoding annotations defined in [selected_annotations]. 
    # for example, if we only have Brain K27ac in the model, then we will flag genes that don't have any non-coding bases covered by this histone marks.
    # And the non-coding bayes factor for these genes would be always 1. Because the null model and the alternative model would be exactly the same. 
    no_feature_flag <- function(data_partition_element){
      no_feature_flag_level2 <-function(data_partition_element_level2){
        sum(data_partition_element_level2[[1]][selected_annotations]^2)
      }
      as.numeric(sum(sapply(data_partition_element, no_feature_flag_level2)) >0) # if there are any bases that have any features that are considered to have significant rr in the model
    }
    
    noncoding_rr_feature_flag <- sapply(data, no_feature_flag)
    # get the names of genes without bases that have informative rr features (as defined in [selected_annotations])
    genes_without_rr_feature <- names(noncoding_rr_feature_flag[noncoding_rr_feature_flag == 0])
    
    
    logBF_noncoding <-data.table(genename = names(logP_Zg1), logBF_noncoding = logP_Zg1 - logP_Zg0)
    # for each gene sum the logBF_noncoding, I did this because there is possibility that a same gene was partitioned into two chunks before collapsing feature data. 
    temp <- by(logBF_noncoding, logBF_noncoding$genename, function(x) {sum(x$logBF_noncoding)})
    logBF_noncoding_list <- append(logBF_noncoding_list, list(data.table(genename = names(temp), logBF_noncoding = as.vector(temp))))
  }
  
  # Sum the log BF of noncoding features over all the studies
  logBF_noncoding <- data.table(genename = logBF_noncoding_list[[1]]$genename, logBF_noncoding = 0)
  for(i in 1:length(logBF_noncoding_list)){
    logBF_noncoding$logBF_noncoding <- logBF_noncoding$logBF_noncoding + logBF_noncoding_list[[i]]$logBF_noncoding
  }
  
  # generate a gene_BF_table to hold up everything, including bayse factors for coding and non-coding mutations
  gene_BF_table <- logBF_noncoding
  # add in coding_BF info
  gene_BF_table <- gene_BF_table[additional_BF, on = "genename"]
  
  # change NA to 0, NA would be generated if no non-coding annotations/mutations of a gene were recorded in  mut_collapsed_data.
  # However, in mut_collapsed_data, a gene would still be included if it has non-coding bases included even in the situation where non of the bases overlap with any epigenomic marks considered or are marked by any other base-level annotation as deletrious mutations
  # For such a gene, need to consider if we want to call it as novel gene, if the combined FDR becomes a little bit less than a cutoff and its coding FDR is a little greater than this cutoff.
  # get the names of genes that don't have any noncoding bases that are ever used in the model, defined in [noncoding_annotations] in [TADA_A_reading_in_annotations]
  genes_without_rr_feature <- union(genes_without_rr_feature, gene_BF_table[is.na(logBF_noncoding)]$genename)
  gene_BF_table[is.na(logBF_noncoding)]$logBF_noncoding = 0
  # add coding and non-coding logBF
  gene_BF_table$logBF_all <- gene_BF_table$logBF_noncoding + gene_BF_table$logBF_coding
  gene_BF_table[,c("BF_noncoding","BF_coding","BF_all")] = exp(gene_BF_table[,c("logBF_noncoding","logBF_coding","logBF_all")])
  gene_BF_table = gene_BF_table[order(gene_BF_table$BF_coding, decreasing = TRUE),]
  gene_BF_table$FDR_coding = Bayesian.FDR(gene_BF_table$BF_coding, TADA_p0)$FDR
  gene_BF_table = gene_BF_table[order(gene_BF_table$BF_all, decreasing = TRUE),]
  gene_BF_table$FDR_all = Bayesian.FDR(gene_BF_table$BF_all, TADA_p0)$FDR
  
  list(gene_BF_table = gene_BF_table, genes_with_no_epi = genes_without_rr_feature)
}


# read in allele-specific annotation, mutation rates, and mutation locations
TADA_A_read_info_AS <- function(mut_file = "../data/Yuen_NM2015_cases_DNM_CADD_gt10.bed", # the bed file with mutations whose allele specific deleterious scores pass a certain threshold
                                         window_file = "../data/Example_windows.bed",
                                         #mutrate_file = "../data/Example_windows.mutrate",
                                         mutrate_scaling_file = "../data/Example_windows_mutrate_scaling_file_for_Yuen_NM2015_cases_DNM.txt",
                                         sample_size = 162,
                                         gene_prior_file = "../data/Example_gene_prior.txt",
                                         report_proportion = 100/18665,
                                         chunk_partition_num =1,
                                         node_n = 6,
                                         AS_mutrate_file_prefix = "../data/nf_and_promoter_whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate", # prefix of CADD mutrate chr by chr. only has regions that are included in the model (e.g., only for non-coding regions within 10kb of TSS)
                                         gene_assignment_mode = "CADD"){
  # [mut_file] is the the bed file with mutations whose allele specific deleterious scores pass a certain threshold
  # [window_file] is the file with genomic windows. Each line represents one window. The columns are "chr", "start", "end" followed by features that might affect local background mutation rate.
  # [sample_size] is the number of individuals
  # [mutrate_scaling_file] is the name of the file that has the scaling factor for each genomic interval in [window_file]. 1st column is site_index, 2nd column is scaling factor of mutation rate. 
  # each mutrate_scaling_file matches to one pair of window_file and mut_file.
  # [gene_prior_file], a file that has prior (derived from posterior and prior)for a gene as a risk gene. 
  # [report_proportion] Choose the top X% TADA genes to estimate RR. 
  # [node_n] is the number of nodes used to run a certain chunk of the code, default is 6
  # [AS_mutrate_file_prefix], is the prefix of allele-specific mutation rate at base-level. e.g., "../data/nf_and_promoter_whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate"
  # mutrate file has four columns, chr, strat, end, and base-level mutrate if [gene_assignment_mode] is "CADD"
  # otherwise has five columns, the extra column is gene assignment. 
  # The mutrate in allele-specific model has four columns, chr, start, end and base-level allele-specific mutation rate. 
  # [estimate_RR] if RR estimation is true. 
  # [gene_assignment_mode] could be "CADD" or "Splicing". if "CADD", a base is assigned to the nearset gene of the genomic window that has the base (Thus, [window_file] has information that determines the assignment).
  # if "splicing", a base is assigned to a gene based on info in the mutrate file specified by [AS_mutrate_file_prefix] 
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  mutrate_scaling <- fread(mutrate_scaling_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- windows[mutrate_scaling, on = "site_index"]
  
  # get the piror probability of genes.
  gene_prior <- fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # merge gene prior info
  coverage_1 <- coverage[gene_prior, on = "genename"]
  coverage_1 <-coverage_1[complete.cases(coverage_1)]
  rm(coverage)
  # select genes based on TADA prior probability and [report_proportion]
  if(report_proportion !=1){
    genes_for_report <- gene_prior[order(gene_prior[,2],decreasing = TRUE),1]
    genes_for_report <- genes_for_report[1:floor(nrow(genes_for_report)*report_proportion)]
    coverage_1 <- coverage_1[genes_for_report, on = "genename"]
    coverage_1 <-coverage_1[complete.cases(coverage_1)]
  }else{
    genes_for_report  <- gene_prior[,1] # choose all the genes in TADA coding table 
  }
  
  # now need to extropolate window mutation file to base level mutation file
  coverage_1 <-coverage_1[,c("chr","start","end","scaling_factor","genename"),with = FALSE]
  
  coverage_1 <- split(coverage_1, coverage_1$chr)
  
  # now get the coordinates of de novo mutations that pass the CADD cutoff threshold
  mut <- fread(mut_file, header = FALSE, sep = "\t")
  
  AS_data_partition <-list()
  
  for(i in 1:length(coverage_1)){
    #write to a bed file, didn't remove mutations that have other annotations for now. 
    fwrite(coverage_1[[i]][,c("chr","start","end","scaling_factor","genename")], paste(prefix, "temp_for_AS_rate.bed",sep = "_"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    AS_mutrate_file <- paste(AS_mutrate_file_prefix,coverage_1[[i]][1,c("chr")], sep = ".")
    command <- paste("bedtools intersect -a ", AS_mutrate_file, " -b ", paste(prefix, "temp_for_AS_rate.bed",sep = "_"), " -wa -wb > ", paste(prefix, "temp_scaling_factor_AS_rate.bed",sep = "_"), sep = "")
    system(command)
    
    if(file.info(paste(prefix, "temp_scaling_factor_AS_rate.bed",sep = "_"))$size == 0) next # this happens when report gene number is small and we are looking at Y chromosome.
    
    mutrate_scaling_AS_rate <- fread(paste(prefix, "temp_scaling_factor_AS_rate.bed",sep = "_"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    # the column names assignment here is different from spidex.
    if(gene_assignment_mode == "CADD"){
      mutrate_scaling_AS_rate <- mutrate_scaling_AS_rate[,c(1,2,3,4,9,8)]
    }else if(gene_assignment_mode == "Splicing"){
      mutrate_scaling_AS_rate <- mutrate_scaling_AS_rate[,c(1,2,3,4,5,9)]
    }
    colnames(mutrate_scaling_AS_rate) <- c("chr","start","end","mutrate","genename","scaling_factor")
    mutrate_scaling_AS_rate$adjusted_mutrate <- mutrate_scaling_AS_rate$mutrate * mutrate_scaling_AS_rate$scaling_factor*2*sample_size # remember to add sample size in. 
    mutrate_scaling_AS_rate <- mutrate_scaling_AS_rate[,sum(.SD$adjusted_mutrate),by = c("chr","start","end","genename"), .SDcols = c("adjusted_mutrate")]
    colnames(mutrate_scaling_AS_rate)[5] <- "adjusted_base_mutrate" # this is the splicing mutation rate summing up allele-specific mutation rates at each base
    
    # get (allele-specific)AS mutation coverage at the base level
    fwrite(mutrate_scaling_AS_rate, paste(prefix,"_temp_AS_for_getting_mutation_coverage.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    command <- paste("bedtools coverage -a ", mut_file, " -b ", paste(prefix,"_temp_AS_for_getting_mutation_coverage.bed", sep = ""), " > ", paste(prefix,"_temp_AS_mut_coverage.bed", sep = ""), sep = "")
    system(command)
    mutrate_scaling_AS_rate <- fread(paste(prefix,"_temp_AS_mut_coverage.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    mutrate_scaling_AS_rate <- mutrate_scaling_AS_rate[,c(1,2,3,4,5,6)]
    colnames(mutrate_scaling_AS_rate) <- c("chr","start","end","genename","adjusted_base_mutrate","mut_count")
    
    mutrate_scaling_AS_rate$feature = 1 # add a feature to everybase, because here, every base does have at least one CADD mutation if mutated
    
    if(nrow(mutrate_scaling_AS_rate) == 0) next # if no genes (in genes for report) on this current chromosome have splicing mutations, exit the current iteration without adding anything to [splicing_data_partition] and go to the next one
    
    partition_AS_feature <- function(pbg){
      # input is one element of the list of partition_by_gene
      pbg_split <- split(pbg, pbg$feature,drop = TRUE)
      feature_combination_number <- length(pbg_split)
      # this function below is different from the function used in dealing with dataset without reading by chunk. Here, prior is not incoporated at this step.
      info_for_each_feature <- function(feature_set){
        list(feature_vector = c(1), sum_mut_rate_count = sum(feature_set$mut_count*log(feature_set$adjusted_base_mutrate)), sum_mut_rate = sum(feature_set$adjusted_base_mutrate), sum_mut_count = sum(feature_set$mut_count), log_fcount = sum(log(factorial(feature_set$mut_count))))
      }
      sapply(pbg_split, info_for_each_feature,simplify = FALSE)
    }
    
    #then partition by gene for the current chromosome
    mutrate_scaling_AS_rate <- split(mutrate_scaling_AS_rate, mutrate_scaling_AS_rate$genename)
    # then partition by feature configuration for each gene in the current chunk
    mutrate_scaling_AS_rate <- sapply(mutrate_scaling_AS_rate, partition_AS_feature, simplify = FALSE)
    # add into [CADD_data_partition]
    AS_data_partition <- append(AS_data_partition, mutrate_scaling_AS_rate)
  }
  system(paste("rm ", prefix, "_temp*", sep = ""))
  
  return(list(AS_info = AS_data_partition))
}


#estimate RR from multiple studies 
TADA_A_RR_estimate_AS <-function(data_list, gene_prior_file, optimization_iteration = 2000){
  #[data_list] is a list of objects [base_info] returned from [TADA_A_reading_in_annotations], each element of the list has compact base-level information for each gene.
  # If you only want to estiamte RR from one single study, just assign data_list = list(one_study$base_info), where one_study is returned from [TADA_A_reading_in_annotations].
  #[gene_prior_file], #[gene_prior_file], a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  # get the piror probability of genes.
  gene_prior <- fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) <- c("genename", "prior")
  gene_prior$genename <- as.character(gene_prior$genename)
  gene_prior$prior <- as.numeric(gene_prior$prior)
  # notice the fr function is different from the function that deals with dataset without partition
  fr <- function(x){
    # now deal with splicing mutations
    AS_rr = x
    AS_cal_logP_Zg1 <- function(data_partition_element){
      AS_cal_logP_Zg1_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[2]]+data_partition_element_level2[[4]]*data_partition_element_level2[[1]]%*%AS_rr-data_partition_element_level2[[3]]*exp(data_partition_element_level2[[1]]%*%AS_rr)-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, AS_cal_logP_Zg1_level2))
    }
    
    AS_cal_logP_Zg0 <- function(data_partition_element){
      AS_cal_logP_Zg0_level2 <-function(data_partition_element_level2){
        data_partition_element_level2[[2]]-data_partition_element_level2[[3]]-data_partition_element_level2[[5]]
      }
      sum(sapply(data_partition_element, AS_cal_logP_Zg0_level2))
    }
    ll_sum1 <- 0
    for(data in data_list){
      AS_logP_Zg1 = sapply(data, AS_cal_logP_Zg1)
      AS_logP_Zg0 = sapply(data, AS_cal_logP_Zg0)
      
      AS_logP_table<-data.table(logP_Zg1 = AS_logP_Zg1, logP_Zg0 = AS_logP_Zg0, genename = names(AS_logP_Zg1))
      
      AS_logP_table <- AS_logP_table[gene_prior, on = "genename"]
      AS_logP_table <- AS_logP_table[complete.cases(AS_logP_table)]
      # get sum of ll of splicing mutations
      sum(by(AS_logP_table, AS_logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
      # sum of ll from all studies
      ll_sum1 <- ll_sum1 + sum(by(AS_logP_table, AS_logP_table$genename, function(x){log(exp(log(x[1,]$prior)+sum(x$logP_Zg1))+exp(log(1-x[1,]$prior)+sum(x$logP_Zg0)))}))
    }
    ll_sum1
  }
  mle <- optim(0.1, fr, method = "Brent", lower = -5, upper = 10, control=list("fnscale"= -1, "maxit" = optimization_iteration), hessian = TRUE)
  
  # get confidence intervals of RR estimates
  fisher_info <- 1/(-diag(mle$hessian))
  prop_sigma <- sqrt(fisher_info)
  rr_estimate <- c(mle$par)
  upper<-rr_estimate+1.96*prop_sigma
  lower<-rr_estimate-1.96*prop_sigma
  rr_report <-data.frame(relative_risk = rr_estimate, lower_bound = lower, upper_bournd = upper)
  
  list(mle = mle, rr_report = rr_report)
}

