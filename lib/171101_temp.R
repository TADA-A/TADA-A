









mut_files = c("../data/Simons_519_new_version_controls_with_allele_info.txt")
window_file = "../data/Example_windows_with_div_score.bed"
mutrate_scaling_files = c("../data/Example_windows_mutrate_with_div_score_scaling_file_for_Simons_519_new_control_DNM.txt")
sample_sizes = c(519)
gene_prior_file = "../data/Example_gene_prior_shuffled.txt"
nonAS_noncoding_annotations = NA
AS_noncoding_annotations = list(c("../other_annotations/coding/171029_synonymous_SNV_altA.bed", "../other_annotations/coding/171029_synonymous_SNV_altC.bed", "../other_annotations/coding/171029_synonymous_SNV_altG.bed", "../other_annotations/coding/171029_synonymous_SNV_altT.bed"))
report_proportion = 50/18665
chunk_partition_num =1
node_n = 6
mutrate_ref_files = c("../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw",
                      "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw",
                      "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw",
                      "../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw")


test1 <- coverage_noncoding_for_base_mutrate_temp[complete.cases(coverage_noncoding_for_base_mutrate_temp)]
test1 <- split(test1, test1$genename)
test1 <- sapply(test1, partition_feature, simplify = FALSE)

data = test1
TADA_A_RR_estimate(data = data, selected_annotations = c(1), gene_prior_file = "../data/Example_gene_prior.txt", optimization_iteration = 2000)

test2 <- coverage_noncoding_for_base_mutrate_temp[complete.cases(coverage_noncoding_for_base_mutrate_temp)]
test2 <- test2[Anno_1_A == 1]
test2 <- split(test2, test2$genename)
test2 <- sapply(test2, partition_feature, simplify = FALSE)

data = test2
TADA_A_RR_estimate(data = data, selected_annotations = c(1), gene_prior_file = "../data/Example_gene_prior.txt", optimization_iteration = 2000)

# all mutations 519 control Alt_A
sum adjusted mutation rates is 79.3  total_mutcount is 177
Seems mutation adjustment is not enough...
sum adjusted synonymous mutation rates is 3.13
synonymous mutation count is 2

$mle
$mle$par
[1] -0.110152

$mle$value
[1] -74.17707

$mle$counts
function gradient 
NA       NA 

$mle$convergence
[1] 0

