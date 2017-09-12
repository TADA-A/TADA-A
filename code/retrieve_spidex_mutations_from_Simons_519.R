library(datat.table)
mutation <- fread("/media/yuwen/F/ASD/160616_Simons_data/cases_with_ref_and_alt.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
AS_ref <- fread("/media/yuwen/F/TADA-A/data/nf_and_promoter_whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr1", header = FALSE, sep = "\t", stringsAsFactors = FALSE)


load("/media/yuwen/F/ASD/data/0703_region_list_080216_data_matrix.Rdata")

mutation_in_data_matrix <- as.data.table(data_matrix[,1:15])
mutation_in_data_matrix <- mutation_in_data_matrix[Prediction == "dnv_proband" & Type == "SNV"]                             
mutation_in_data_matrix$Start <- mutation_in_data_matrix$Start - 1
fwrite(mutation_in_data_matrix[,1:3], "/media/yuwen/F/TADA-A/data/Simons_519_new_version_cases.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
mutation_in_data_matrix_cadd_gt_10 <- mutation_in_data_matrix[CADD13_PHRED > 10 & !is.na(CADD13_PHRED)]
fwrite(mutation_in_data_matrix_cadd_gt_10[,1:3], "/media/yuwen/F/TADA-A/data/Simons_519_new_version_cases_CADD_gt10.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

temp <- data.table(mutation_in_data_matrix[,1], mutation_in_data_matrix[,3], index = seq(1:dim(mutation_in_data_matrix)[1]), mutation_in_data_matrix[,4:5])
fwrite(temp, "/media/yuwen/F/TADA-A/data/Simons_519_new_version_cases_for_getting_spidex.txt", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

temp1 <- fread("/media/yuwen/F/TADA-A/data/Simons_519_new_version_cases_for_getting_spidex.results", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
temp1 <- temp1[!is.na(dpsi_zscore) & dpsi_zscore < -1.416]

test <- merge(temp1, temp, by = "index")

test<- data.table(test$Chrom, test$End -1 , test$End)
fwrite(test, "/media/yuwen/F/TADA-A/data/Simons_519_new_version_cases_spidex_lower_10pct.bed", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
