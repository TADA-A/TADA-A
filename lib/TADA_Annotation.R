

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
  command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", mut_file, " -b ", paste(prefix, "_temp_windows.bed", sep = ""), " > ", paste(prefix,"_temp_coverage.bed", sep = ""), sep = "")
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
                                                   "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw"),
                             MPI = 1){
  
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
  # [MPI] is the index that will add to temp files, useful when running multipe processes at one time
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, MPI, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  #mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- windows
  # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  for(i in 1:length(mut_files)){
    mutrate_scaling <- fread(mutrate_scaling_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    system(paste("echo \"Finished reading mutrate scaling file ", mutrate_scaling_files[i], ".\"", sep = ""))
    system("date")
    coverage <- coverage[mutrate_scaling, on = "site_index"]
    coverage <- coverage[complete.cases(coverage)] # release memory
    colnames(coverage)[length(colnames(coverage))] <- paste("scaling_factor_study_", i, sep = "")
    rm(mutrate_scaling) # release memory
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
    data.frame(table_row[1], start, start+1, paste(table_row["chr"],table_row["genename"],start,sep = "_"), table_row["ID"])
  }
  
  options(warn=-1)
  if(node_n != 1){
    cl <- makeCluster(node_n)
  }
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  
  # get nonAS feature number
  if(is.na(nonAS_noncoding_annotations[1])){
    nonAS_feature_number <- 0 
  }else{
    nonAS_feature_number <- length(nonAS_noncoding_annotations)
  }
  # get AS feature number
  if(is.na(AS_noncoding_annotations[1])){
    AS_feature_number <- 0 
  }else{
    AS_feature_number <- length(AS_noncoding_annotations)
  }
  
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
    if(node_n != 1){
      coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage[[i]], 1, window_expansion))
    }else{
      coverage_noncoding_for_base_mutrate <- rbindlist(apply(coverage[[i]], 1, window_expansion))
    }
    system(paste("echo \"Finished partitioning base-level coordinates data at Round ", i, ".\"", sep = ""))
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
      system(paste("echo \"Finished obtaining base-level mutation rate for alt allele ", alt_letters[j], ".\"", sep = ""))
      system("date")
    }
    
    # read in non allele-specific epigenomic annotations
    epi_ID = 1
    if (!is.na(nonAS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in nonAS_noncoding_annotations){
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, sep = "_"))
        coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        system(paste("echo \"Finished reading non-allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
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
          command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi[k], " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
          system(command)
          base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
          base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
          colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, alt_letters[k], sep = "_"))
          coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        }
        system(paste("echo \"Finished reading allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
        system("date")
        epi_ID <- epi_ID + 1
      }
    }
    
    
    # now for each study, read in data, and collapse data based on noncoding annotation configuration
    for(m in 1:length(mut_files)){
      mut <- fread(mut_files[m], header = FALSE, sep = "\t", stringsAsFactors = FALSE)
      for(letter in alt_letters){
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }
        # a very important step here is removing bases that have adjusted_mutrate_base 0. This happens when the allele of the mutrate we are using is just the reference allele. 
        # By doing this, we make the computation of likelihood valid, also we automatically removed bases with nonAS annotations but with mutant allele the same with ref allele at this step
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage_noncoding_for_base_mutrate_temp[[paste("base_mutrate_alt_", letter, sep = "")]] !=0]
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage[[i]][,c(paste("scaling_factor_study_", m, sep = ""), "genename","ID"), with = FALSE],on = "ID"]
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- coverage_noncoding_for_base_mutrate_temp[, paste("base_mutrate_alt_", letter, sep = ""), with = FALSE] * coverage_noncoding_for_base_mutrate_temp[, paste("scaling_factor_study_", m, sep = ""), with = FALSE] * 2 * sample_sizes[m]
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(1, nonAS_feature_number), sep = "")), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[, c("base_ID", "adjusted_base_mutrate", "genename", paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = "")), with = FALSE]
        }
        mut_allele <- mut[V5 == letter]
        # the subsequent two lines of code are used to prevent a outputing bug when using fwrite. (85000000 to be written as 8.5e7)
        mut_allele$V2 <- as.integer(mut_allele$V2)
        mut_allele$V3 <- as.integer(mut_allele$V3)
        fwrite(mut_allele, paste(prefix, "_temp_mut_allele.bed", sep = ""), col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", paste(prefix, "_temp_mut_allele.bed", sep = ""), " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""),sep = "")
        system(command)
        base_with_mut <- fread(paste(prefix,"_temp_for_mutrate_overlap_mut_allele.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_with_mut <- base_with_mut[,c("V4","V5"), with = FALSE]
        colnames(base_with_mut) <- c("base_ID", "mut_count")
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[base_with_mut, on = "base_ID"]
        
        # have to collpase data at this point, otherwise, the I will run out of RAM if process 1000 genes every time. 
        anno_count <- rep(0, nrow(coverage_noncoding_for_base_mutrate_temp))
        for(p in 1:feature_number){
          anno_count <- anno_count + coverage_noncoding_for_base_mutrate_temp[[(4 + p -1)]]
        }
        # remove rows(bases) that don't have any non-coding features, this could save a lot of RAM, so I could use smaller partition number which would greatly accelerate speed when read in data for all genes.
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[anno_count >0,]
        
        # first partition by gene for the current chunk
        coverage_noncoding_for_base_mutrate_temp<- split(coverage_noncoding_for_base_mutrate_temp, coverage_noncoding_for_base_mutrate_temp$genename)
        # then partition by feature configuration for each gene in the current chunk
        coverage_noncoding_for_base_mutrate_temp <- sapply(coverage_noncoding_for_base_mutrate_temp, partition_feature, simplify = FALSE)
        # add compact data
        data_partition <- append(data_partition, coverage_noncoding_for_base_mutrate_temp)
        rm(coverage_noncoding_for_base_mutrate_temp) # release memory
        system(paste("echo \"Finished read in mutation data and make them into the compact format for Study ", m, " and allele ", letter, ".\"", sep = ""))
        system("date")
      }
    }
  }
  
  if(node_n != 1){
    stopCluster(cl)
  }
  options(warn = 0)
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste("echo \"Temp files cleaned and data recording finished!\""))
  system("date")
  return(list(base_info = data_partition))
}


TADA_A_RR_estimate <-function(data, selected_annotations, gene_prior_file, optimization_iteration = 2000, mode = "regular"){
  #[data] is the [base_info] returned from [TADA_A_reading_in_annotations], which contains all the allele specific data across all studies
  #[selected_annotations] is a vector indicating non-coding annotations whose RRs need to be estimated. e.g., c(2,3) means that the 2nd and 3rd annotations in the [noncoding_annotations] argument of [TADA_A_reading_in_annotations] will have their RRs estimated.
  #[gene_prior_file], #[gene_prior_file], a file that has the prior probability of each gene as a risk gene. e.g., "../data/Example_gene_prior.txt".
  #[optimization_iteration] is the number of iterations that optim() will perform to estimate RRs.
  #[mode] is "regular", or "single_fast". "single_fast" is used when estimating RR from only one annotation ([data] only recoreded one annotation) of lots of genes (e.g., all genes), would be 5 times faster.
  
  # get the piror probability of genes.
  gene_prior = fread(gene_prior_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  colnames(gene_prior) = c("genename", "prior")
  gene_prior$prior <- as.numeric(gene_prior$prior)
  gene_prior$genename = as.character(gene_prior$genename)
  if(mode == "single_fast"){
    further_collapsed_data <- list()
    genes_with_annotation <- unique(names(data))
    
    for(i in 1:length(genes_with_annotation)){
      further_collapsed_data <- append(further_collapsed_data, list(list(`1` = list(feature_vector = 1 , sum_mut_rate_count= 0, sum_mut_rate = 0, sum_mut_count = 0, log_fcount = 0))))
    }
    names(further_collapsed_data) <- unique(names(data))
    
    for(i in 1:length(data)){
      Part_of_data = data[i]
      further_collapsed_data[[names(Part_of_data)]][[1]]$sum_mut_rate_count <- further_collapsed_data[[names(Part_of_data)]][[1]]$sum_mut_rate_count + Part_of_data[[1]][[1]]$sum_mut_rate_count
      further_collapsed_data[[names(Part_of_data)]][[1]]$sum_mut_rate <- further_collapsed_data[[names(Part_of_data)]][[1]]$sum_mut_rate + Part_of_data[[1]][[1]]$sum_mut_rate
      further_collapsed_data[[names(Part_of_data)]][[1]]$sum_mut_count <- further_collapsed_data[[names(Part_of_data)]][[1]]$sum_mut_count + Part_of_data[[1]][[1]]$sum_mut_count
      further_collapsed_data[[names(Part_of_data)]][[1]]$log_fcount <- further_collapsed_data[[names(Part_of_data)]][[1]]$log_fcount + Part_of_data[[1]][[1]]$log_fcount
    }
    data <- further_collapsed_data
  }
  
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


# Bayesian FDR control (PMID:19822692, Section2.3)

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


# This function is for calculating Bayes factors for based on noncoding annotations and then combine with coding Bayes factors in order to predict risk genes.
# calculate the Bayes factor for each gene based on informative noncodng annotations. 
TADA_A_get_BFs <- function(data,selected_annotations, rr, additional_BF_file, TADA_p0 = 0.94){
  # [data] is the [base_info] returned from [TADA_A_read_info], each element of the list has compact base-level information for each gene. 
  # [selected_annotations] determines which features in the [data] are going to be used for calculating bayes factors. e.g., if [selected_annotations] is c(1,3). Then the first and third element of feature_vector in the [base_info] of [data] will be selected.
  # [rr] is the estimated relative risks of [selected_annotations]. RRs are estimated from [TADA_A_RR_estimate]. 
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
  genes_with_rr_feature <- unique(names(noncoding_rr_feature_flag[noncoding_rr_feature_flag != 0])) # a gene may present in multiple elements in [data], so use unique to get a unique list of genes.
  genes_without_rr_feature <- additional_BF[!is.element(additional_BF$genename, genes_with_rr_feature),]$genename  
  
  logBF_noncoding <-data.table(genename = names(logP_Zg1), logBF_noncoding = logP_Zg1 - logP_Zg0)
  # for each gene sum the logBF_noncoding, I did this because there is possibility that a same gene was partitioned into two chunks before collapsing feature data. 
  logBF_noncoding <- by(logBF_noncoding, logBF_noncoding$genename, function(x) {sum(x$logBF_noncoding)})
  logBF_noncoding <- data.table(genename = names(logBF_noncoding), logBF_noncoding = as.vector(logBF_noncoding))
  
  # generate a gene_BF_table to hold up everything, including bayse factors for coding and non-coding mutations
  gene_BF_table <- logBF_noncoding
  # add in coding_BF info
  gene_BF_table <- gene_BF_table[additional_BF, on = "genename"]
  
  
  # change NA to 0, NA would be generated if no non-coding annotations/mutations of a gene were recorded in  mut_collapsed_data.
  # However, in mut_collapsed_data, a gene would still be included if it has non-coding bases included even in the situation where non of the bases overlap with any epigenomic marks considered or are marked by any other base-level annotation as deletrious mutations
  # For such a gene, need to consider if we want to call it as novel gene, if the combined FDR becomes a little bit less than a cutoff and its coding FDR is a little greater than this cutoff.
  # get the names of genes that don't have any noncoding bases that are never used in the model, defined in [noncoding_annotations] in [TADA_A_read_info]
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



TADA_A_DNM_generator <- function(window_file = "../data/Example_windows.bed",
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
                                                       "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw"),
                                 rr = c(1, 0, 0.5),
                                 output_allele_info_files,
                                 output_bed_files,
                                 output_risk_genes_file = "temp.riskgenes",
                                 compact_mut_output = NA,
                                 MPI = 1){
  
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
  # [rr] is the log-relative risk of DNMs.
  # [output_allele_info_files] is a vector of files containing coordinates of mutations and ref/alt alleles
  # [output_bed_files] is a vector of files containing coordinates only of mutations.
  # [output_risk_genes_file] is a string indicating the output file containing the names of risk genes in a column.
  # [compact_mut_output] if set to be True, a RDS with mutation data in a compact form will be generated. Useful for estimating RR and claculating Bayes factor
  # [MPI] is the index that will add to temp files, useful when running multipe processes at one time
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, MPI, sep = "")
  
  # make a tmp folder for tmp files
  system("mkdir -p tmp")
  
  #mut <- fread(mut_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  windows <- fread(window_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  coverage <- windows
  # the number of genomic windows in [mutrate_scaling] is less than the number of windows in [windows] because there are a few windows with mutration rate equal to 0, and thus removed.
  for(i in 1:length(mutrate_scaling_files)){
    mutrate_scaling <- fread(mutrate_scaling_files[i], header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    system(paste("echo \"Finished reading mutrate scaling file ", mutrate_scaling_files[i], ".\"", sep = ""))
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
  genes_for_report_with_prior <- gene_prior[order(gene_prior[,2],decreasing = TRUE),]
  genes_for_report_with_prior <- genes_for_report_with_prior[1:floor(nrow(genes_for_report_with_prior)*report_proportion),]
  # randomly assign risk genes based on prior probability
  genes_for_report_with_prior$risk <- rbinom(nrow(genes_for_report_with_prior), 1, genes_for_report_with_prior$prior)
  risk_genes <- genes_for_report_with_prior[risk == 1,]$genename
  nonrisk_genes <- genes_for_report_with_prior[risk == 0,]$genename
  coverage_risk <- coverage[is.element(genename, risk_genes)]
  coverage_risk$risk <- 1
  coverage_nonrisk <- coverage[is.element(genename, nonrisk_genes)]
  coverage_nonrisk$risk <- 0 
  
  coverage <- rbind(coverage_risk, coverage_nonrisk)
  
  # now need to extropolate window mutation file to base level mutation file
  coverage <-coverage[,c("chr","start","end",paste("scaling_factor_study_", seq(1,length(mutrate_scaling_files)), sep = ""),"genename","risk"),with = FALSE]
  
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
  
  
  bed_extension <- function(table_row){
    data.frame(chr = rep(table_row["chr"], table_row["mut_count"]), start = rep(table_row["start"], table_row["mut_count"]), end = rep(table_row["end"], table_row["mut_count"]))
  }
  
  options(warn=-1)
  if(node_n != 1){
    cl <- makeCluster(node_n)
  }
  # use parallel computing and rbindlist to increase the speed by thousands of times. 
  environment(window_expansion) <- .GlobalEnv
  environment(bed_extension) <- .GlobalEnv
  # get nonAS feature number
  if(is.na(nonAS_noncoding_annotations[1])){
    nonAS_feature_number <- 0 
  }else{
    nonAS_feature_number <- length(nonAS_noncoding_annotations)
  }
  # get AS feature number
  if(is.na(AS_noncoding_annotations[1])){
    AS_feature_number <- 0 
  }else{
    AS_feature_number <- length(AS_noncoding_annotations)
  }
  
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
  
  
  alt_letters <- c("A","C","G","T")
  mut_output_list <- list()
  data_partition <- list()
  for(m in 1:length(mutrate_scaling_files)){
    mut_output_list[[m]] <- data.table(chr = character(), start = integer(), end = integer(), ref = character(), alt = character())
  }
  
  # here only deals with the noncoding parts that are within 10kb of TSSs of genes. Will deal with 
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    if(node_n != 1){
      coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage[[i]], 1, window_expansion))
    }else{
      coverage_noncoding_for_base_mutrate <- rbindlist(apply(coverage[[i]], 1, window_expansion))
    }
    system(paste("echo \"Finished partitioning base-level coordinates data at Round ", i, ".\"", sep = ""))
    system("date")
    colnames(coverage_noncoding_for_base_mutrate) <- c("chr","start","end","base_ID","ID")
    coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[!duplicated(base_ID)]
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
      system(paste("echo \"Finished obtaining base-level mutation rate for alt allele ", alt_letters[j], ".\"", sep = ""))
      system("date")
    }
    
    # read in non allele-specific epigenomic annotations
    epi_ID = 1
    if (!is.na(nonAS_noncoding_annotations)[1]){ # then epigenomic_marks must be a vector of epigenomic bed files that need to be compard with the mutation data
      for(epi in nonAS_noncoding_annotations){
        command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
        system(command)
        base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
        base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
        colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, sep = "_"))
        coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        system(paste("echo \"Finished reading non-allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
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
          command <- paste("../external_tools/bedtools-2.17.0/bin/bedtools coverage -a ", epi[k], " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
          system(command)
          base_in_epi <- fread(paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
          base_in_epi <- base_in_epi[,c("V4","V5"), with = FALSE]
          colnames(base_in_epi) <- c("base_ID", paste("Anno",epi_ID, alt_letters[k], sep = "_"))
          coverage_noncoding_for_base_mutrate <- coverage_noncoding_for_base_mutrate[base_in_epi, on = "base_ID"]
        }
        system(paste("echo \"Finished reading allele specific noncoding annotations ", epi_ID, ".\"", sep = ""))
        system("date")
        epi_ID <- epi_ID + 1
      }
    }
    
    # now for each sample size, generate mut data.
    for(m in 1:length(mutrate_scaling_files)){
      mut_output_temp <- data.table(chr = character(), start = integer(), end = integer(), ref = character(), alt = character())
      for(letter in alt_letters){
        if(!is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = ""), "chr", "start", "end"), with = FALSE]
        }else if(!is.na(nonAS_noncoding_annotations)[1] & is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(1, nonAS_feature_number), sep = ""), "chr", "start", "end"), with = FALSE]
        }else if(is.na(nonAS_noncoding_annotations)[1] & !is.na(AS_noncoding_annotations)[1]){
          coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate[, c("base_ID", "ID",paste("base_mutrate_alt_", letter, sep = ""), paste("Anno_", seq(nonAS_feature_number +1, nonAS_feature_number + AS_feature_number), "_", letter ,sep = ""), "chr", "start", "end"), with = FALSE]
        }
        # a very important step here is removing bases that have adjusted_mutrate_base 0. This happens when the allele of the mutrate we are using is just the reference allele. 
        # By doing this, we make the computation of likelihood valid, also we automatically removed bases with nonAS annotations but with mutant allele the same with ref allele at this step
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage_noncoding_for_base_mutrate_temp[[paste("base_mutrate_alt_", letter, sep = "")]] !=0]
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[coverage[[i]][,c(paste("scaling_factor_study_", m, sep = ""), "genename","ID","risk"), with = FALSE],on = "ID"]
        # removing sites with allele-specific mutrate == 0 will result in incomplete rows after merging.
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[complete.cases(coverage_noncoding_for_base_mutrate_temp)]
        # derive the adjusted mutation rate for bases of risk genes and nonrisk genes.
        # first multiplies sample size and scaling factor, adjusted_base_mutrate will be used for generating data in a compact format.
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- coverage_noncoding_for_base_mutrate_temp[, paste("base_mutrate_alt_", letter, sep = ""), with = FALSE] * coverage_noncoding_for_base_mutrate_temp[, paste("scaling_factor_study_", m, sep = ""), with = FALSE] * 2 * sample_sizes[m]
        # Then multiplies effect size after taking exponential, to generate a mutrate that has taken into account RRs, will be used for simulating mutations
        # 
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate_rr <- exp(as.matrix(coverage_noncoding_for_base_mutrate_temp[,4 : (3 + feature_number)]) %*% rr * coverage_noncoding_for_base_mutrate_temp$risk) * coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate
        # generate random allele-specific mutations based on adjusted mutation rate.
        coverage_noncoding_for_base_mutrate_temp$mut_count <- rpois(nrow(coverage_noncoding_for_base_mutrate_temp), coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate_rr)
        
        if(!is.na(compact_mut_output)){
          # have to collpase data at this point, otherwise, the I will run out of RAM if process 1000 genes every time. 
          anno_count <- rep(0, nrow(coverage_noncoding_for_base_mutrate_temp))
          for(p in 1:feature_number){
            anno_count <- anno_count + coverage_noncoding_for_base_mutrate_temp[[(4 + p -1)]]
          }
          # remove rows(bases) that don't have any non-coding features, this could save a lot of RAM, so I could use smaller partition number which would greatly accelerate speed when read in data for all genes.
          coverage_noncoding_for_base_mutrate_temp1 <- coverage_noncoding_for_base_mutrate_temp[anno_count >0,]
          
          # first partition by gene for the current chunk
          coverage_noncoding_for_base_mutrate_temp1<- split(coverage_noncoding_for_base_mutrate_temp1, coverage_noncoding_for_base_mutrate_temp1$genename)
          # then partition by feature configuration for each gene in the current chunk
          coverage_noncoding_for_base_mutrate_temp1 <- sapply(coverage_noncoding_for_base_mutrate_temp1, partition_feature, simplify = FALSE)
          # add compact data
          data_partition <- append(data_partition, coverage_noncoding_for_base_mutrate_temp1)
          rm(coverage_noncoding_for_base_mutrate_temp1) # release memory
          system(paste("echo \"Finished read in mutation data and make them into the compact format for Study ", m, " and allele ", letter, ".\"", sep = ""))
          system("date")
        }
        # then generate bed file format, this is a pseudobed file as I will put N at the 4th column as the ref allele is not relevant and mutant allele at the last column.
        # remove rows(bases) that don't have any DNMs
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[mut_count >0,]
        if(node_n != 1){
          mut_bed <- rbindlist(parApply(cl, coverage_noncoding_for_base_mutrate_temp, 1, bed_extension))
        }else{
          mut_bed <- rbindlist(apply(coverage_noncoding_for_base_mutrate_temp, 1, bed_extension))
        }
        mut_bed$ref <- "N"
        mut_bed$alt <- letter
        mut_bed$start <- as.integer(as.character(mut_bed$start))
        mut_bed$end <- as.integer(as.character(mut_bed$end))
        mut_output_temp <- rbind(mut_output_temp, mut_bed)
      }
      mut_output_list[[m]] <- rbind(mut_output_list[[m]], mut_output_temp)
    }
  }
  
  if(node_n != 1){
    stopCluster(cl)
  }
  options(warn = 0)
  # remove temporary files
  system(paste("rm ", prefix, "_temp*", sep = ""))
  print(paste("echo \"Temp files cleaned and mutation generation finished finished!\""))
  system("date")
  # now write mutation file and mutation with allele info file
  for(m in 1:length(mutrate_scaling_files)){
    fwrite(mut_output_list[[m]], output_allele_info_files[m], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
    fwrite(mut_output_list[[m]][,1:3], output_bed_files[m], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  }
  write.table(risk_genes, output_risk_genes_file, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  if(!is.na(compact_mut_output)){
    saveRDS(data_partition, compact_mut_output)
  }
}



retrieve_mutation_info <- function(data_partition, genename, noncoding_annotation){
  #[data_partition] is the [base_info] returned from [TADA_A_read_info]
  #[genename] is the genename to search for
  # [noncoding_annotation] is a vector containing the annotation names that stored in data_partition
  
  data_of_a_gene <- data_partition[names(data_partition) == genename]
  for(i in 1:length(data_of_a_gene)){
    for(j in 1:length(data_of_a_gene[[i]])){
      if(data_of_a_gene[[i]][[j]]$sum_mut_count > 0){
        print(noncoding_annotation[data_of_a_gene[[i]][[j]]$feature_vector * seq(1,length(noncoding_annotation))])
        print(data_of_a_gene[[i]][[j]]$sum_mut_count)
      }
    }
  }
}
