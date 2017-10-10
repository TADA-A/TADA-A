# function to get the number of coding mutations from window files that I have used to train mutation rate models. 
get_coding_mut_number<-function(mut_info, ref_alt_allele, window_file, annovar_folder, phenotype, geneset = "no"){
  # [mut_info] is the mutation information in the data.frame {mutation} in Rdata {*transformed_for_old_code.Rdata}
  # [ref_alt_allele] is the information of ref/mutation alleles stored as a data.frame in {ref_alt_allele} in Rdata {*transformed_for_old_code.Rdata}
  # [window_file] is the file with genome partitioned into windows. e.g., ../other_annotation/epigenomic_annotation/Whole_genome.promoter_yanyu_10kb_exons.50bp_window.bed
  # [annovar_folder] is the folder that has annovar.pl. e.g., /media/yuwen/Elements/ANNOVAR/annovar/
  # [geneset] if "no", then mutation counts of different types for all the genes will be returned. Otherwise, could specify strigent_ASD, relaxed_ASD or other genesets. 
  # use Annovar to annotate coding mutations. 
  # [phenotype] is the phenotype that I want to count, for example "ASD" or "control"
  
  prefix=system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline
  temp_for_annovar = data.frame(mut_info[,1:3], ref_alt_allele, mut_info$mut_type, mut_info$index, mut_info$phenotype)
  temp_for_annovar = temp_for_annovar[temp_for_annovar$mut_info.phenotype == phenotype,1:7]
  write.table(temp_for_annovar, paste(prefix, "_for_annovar_temp.txt",sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  
  # only keep mutations that overlap with coding portions from [window_file]
  command = paste("grep coding ", window_file, " > ", paste(prefix, "_coding_window_temp.bed", sep = ""), sep = "")
  system(command)
  command = paste("bedtools intersect -a ",paste(prefix, "_for_annovar_temp.txt",sep = ""), " -b ", paste(prefix, "_coding_window_temp.bed", sep = ""), " -wa > ",paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), sep = "")
  system(command)
  command = paste("perl ", annovar_folder, "annotate_variation.pl -geneanno -buildver hg19 ", paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""), " ", annovar_folder, "humandb/", sep = "")
  system(command)
  coding_anno = read.delim(paste(paste(prefix, "_for_annovar_temp_in_coding.txt",sep = ""),".exonic_variant_function", sep = ""), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  coding_gene = unlist(lapply(strsplit(coding_anno[,3],split = ":"),`[[`,1))
  coding_mut_for_gene = data.frame(mut_type_raw = coding_anno[,2], genename = coding_gene)
  coding_mut_for_gene = data.frame(coding_mut_for_gene, mut_type = "not_defined")
  coding_mut_for_gene$mut_type = as.character( coding_mut_for_gene$mut_type)
  coding_mut_for_gene[coding_mut_for_gene$mut_type_raw == "synonymous SNV",]$mut_type = "synonymous"
  coding_mut_for_gene[coding_mut_for_gene$mut_type_raw == "nonsynonymous SNV",]$mut_type = "nonsynonymous"
  coding_mut_for_gene[coding_mut_for_gene$mut_type_raw == "stopgain" | coding_mut_for_gene$mut_type_raw == "frameshift deletion"| coding_mut_for_gene$mut_type_raw == "frameshift insertion" ,]$mut_type = "LoF"
  
  if(geneset != "no"){
    coding_mut_for_gene = coding_mut_for_gene[is.element(coding_mut_for_gene$gene, geneset), ]
  }
  #exon_mut_type = rep("unknown", length(mutation[,1]))
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "frameshift deletion",10])] = "frameshift deletion"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "frameshift insertion",10])] = "frameshift insertion"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "nonframeshift deletion",10])] = "nonframeshift deletion"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "nonsynonymous SNV",10])] = "nonsynonymous"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "stopgain",10])] = "stopgain"
  #exon_mut_type[is.element(mutation$index, coding_anno[coding_anno[,2] == "synonymous SNV",10])] = "synonymous"
  
  system(paste("rm ", prefix, "_for_annovar_temp.txt*",sep = ""))
  
  mutation_number = list(syn_number = nrow(coding_mut_for_gene[coding_mut_for_gene[,3] == "synonymous",]),
                         nonsyn_number = nrow(coding_mut_for_gene[coding_mut_for_gene[,3] == "nonsynonymous",]),
                         lof_number = nrow(coding_mut_for_gene[coding_mut_for_gene[,3] == "LoF",]))
  list(mutation_number = mutation_number, coding_mut_for_gene = coding_mut_for_gene)
}
