# This script is used to generate allele-specific base level mutation rate

# Step 1: prepare extended genomic windows file
sed '1d' ../data/Example_windows.bed | awk {'print $1"\t"$2-1"\t"$3+1"\t"$4'} > ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed


# Step 2: Get the nucleotide sequence of each interval, in tab format
bedtools getfasta -fi ../other_annotations/genome_build/hg19.fasta -bed ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed -fo ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta -tab

# Step 3: Use the output file to extract tri-nuleotide sequence of each base within the window intervals.
python ../external_tools/tri_extract_for_TADA-A.py ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta > ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri
rm ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta
# Step 4: Use .tri file as an input file to get the allele-specific mutation rate. Asssing mutation rate 
python ../external_tools/MutRateBase_for_TADA-A_v2.py ../other_annotations/Mark_Daly_mutrate/fordist_1KG_mutation_rate_table.txt ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri > ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.mutrate
rm ../data/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri

# Step 5: generte mutrate file base on the alternative nucleotide, to be filled.


# Step 6: Generate alternative-allele-specific Wig file

sh ../external_tools/base_mutarate_to_wiggle_file.sh ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate
sh ../external_tools/base_mutarate_to_wiggle_file.sh ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate
sh ../external_tools/base_mutarate_to_wiggle_file.sh ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate
sh ../external_tools/base_mutarate_to_wiggle_file.sh ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate

# Step 7: tranform Wig file to bigwig file
# 
../external_tools/wigToBigWig ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.wiggle ../other_annotations/genome_build/hg19.genome ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_A.mutrate.bw
../external_tools/wigToBigWig ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.wiggle ../other_annotations/genome_build/hg19.genome ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_T.mutrate.bw
../external_tools/wigToBigWig ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.wiggle ../other_annotations/genome_build/hg19.genome ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_C.mutrate.bw
../external_tools/wigToBigWig ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.wiggle ../other_annotations/genome_build/hg19.genome ../other_annotations/Mark_Daly_mutrate/Example_windows_extended_1bp_for_getting_base_level_mutrate.bed.fasta.tri.alt_G.mutrate.bw