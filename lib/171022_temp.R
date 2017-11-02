window_file = "../data/Example_windows.bed"
mutrate_scaling_files = c("../data/Example_windows_mutrate_scaling_file_for_Yuen_NM2015_cases_DNM.txt","../data/Example_windows_mutrate_scaling_file_for_Kong_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Wu_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Jiang_cases_DNM.txt", "../data/Example_windows_mutrate_scaling_file_for_Michaelson_cases_DNM.txt")
sample_sizes = c(1620, 780, 320, 320, 100)
gene_prior_file = "../data/Example_gene_prior.txt"
nonAS_noncoding_annotations = c("../data/Noonan_brain_roadmap_union_within_10kb_and_promoter_no_utr.bed","../data/Epigenome_E081_E082_intersection__within_10kb_and_promoter_no_utr.bed","../data/Encode_DHS_union_within_10kb_and_promoter_no_utr.bed")
AS_noncoding_annotations = NA
report_proportion = 100/18665
chunk_partition_num = 1
node_n = 6
mutrate_ref_files = c("../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                      "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                      "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                      "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw")
rr = c(1, 0, 0.5)













# This function randomly generate de novo mutations

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
                                 output_bed_files){
  
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
  # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  prefix <- paste("tmp/", prefix, sep = "")
  
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
  cl <- makeCluster(node_n)
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
  for(m in 1:length(mutrate_scaling_files)){
    mut_output_list[[m]] <- data.table(chr = character(), start = integer(), end = integer(), ref = character(), alt = character())
  }
  
  # here only deals with the noncoding parts that are within 10kb of TSSs of genes. Will deal with 
  for(i in 1:chunk_partition_num){
    # split coverage_noncoding to 10 parts, and for each part, generate feature table (which will be used for )
    coverage_noncoding_for_base_mutrate <- rbindlist(parApply(cl, coverage[[i]], 1, window_expansion))
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
        command <- paste("bedtools coverage -a ", epi, " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
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
          command <- paste("bedtools coverage -a ", epi[k], " -b ", paste(prefix, "_temp_for_mutrate.bed", sep = ""),  " > ", paste(prefix,"_temp_for_mutrate_overlap_epi.bed", sep = ""),sep = "")
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
        # first multiplies sample size and scaling factor
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- coverage_noncoding_for_base_mutrate_temp[, paste("base_mutrate_alt_", letter, sep = ""), with = FALSE] * coverage_noncoding_for_base_mutrate_temp[, paste("scaling_factor_study_", m, sep = ""), with = FALSE] * 2 * sample_sizes[m]
        # Then multiplies effect size after taking exponential 
        coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate <- exp(as.matrix(coverage_noncoding_for_base_mutrate_temp[,4 : (3 + feature_number)]) %*% rr * coverage_noncoding_for_base_mutrate_temp$risk) * coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate
        # generate random allele-specific mutations based on adjusted mutation rate.
        coverage_noncoding_for_base_mutrate_temp$mut_count <- rpois(nrow(coverage_noncoding_for_base_mutrate_temp), coverage_noncoding_for_base_mutrate_temp$adjusted_base_mutrate)
        # then generate bed file format, this is a pseudobed file as I will put N at the 4th column as the ref allele is not relevant and mutant allele at the last column.
        # remove rows(bases) that don't have any DNMs
        coverage_noncoding_for_base_mutrate_temp <- coverage_noncoding_for_base_mutrate_temp[mut_count >0,]
        mut_bed <- rbindlist(parApply(cl, coverage_noncoding_for_base_mutrate_temp, 1, bed_extension))
        mut_bed$ref <- "N"
        mut_bed$alt <- letter
        mut_bed$start <- as.integer(as.character(mut_bed$start))
        mut_bed$end <- as.integer(as.character(mut_bed$end))
        mut_output_temp <- rbind(mut_output_temp, mut_bed)
      }
      mut_output_list[[m]] <- rbind(mut_output_list[[m]], mut_output_temp)
    }
  }
  
  stopCluster(cl)
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
}
