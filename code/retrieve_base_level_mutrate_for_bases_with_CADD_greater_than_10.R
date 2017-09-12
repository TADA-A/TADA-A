bedtools getfasta -fi ../genome_build/hg19.fasta -bed whole_genome_SNVs.gt10.merged.bed -fo whole_genome_SNVs.gt10.merged.bed.fasta
python ../Yi_allele_specific_mutation_rate/Mark_Daly_mutrate/tri_extract_changed_by_YL.py whole_genome_SNVs.gt10.merged.bed.fasta > whole_genome_SNVs.gt10.merged.bed.fasta.tri

python ../Yi_allele_specific_mutation_rate/Mark_Daly_mutrate/MutRateBase_v2.py ../Yi_allele_specific_mutation_rate/Mark_Daly_mutrate/fordist_1KG_mutation_rate_table.txt whole_genome_SNVs.gt10.merged.bed.fasta.tri > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate 

grep -P "chr1\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr1
grep -P "chr2\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr2
grep -P "chr3\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr3
grep -P "chr4\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr4
grep -P "chr5\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr5
grep -P "chr6\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr6
grep -P "chr7\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr7
grep -P "chr8\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr8
grep -P "chr9\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr9
grep -P "chr10\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr10
grep -P "chr11\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr11
grep -P "chr12\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr12
grep -P "chr13\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr13
grep -P "chr14\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr14
grep -P "chr15\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr15
grep -P "chr16\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr16
grep -P "chr17\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr17
grep -P "chr18\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr18
grep -P "chr19\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr19
grep -P "chr20\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr20
grep -P "chr21\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr21
grep -P "chr22\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr22
grep -P "chrX\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chrX
grep -P "chrY\t" whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate > whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chrY

for i in whole_genome_SNVs.gt10.merged.bed.fasta.tri.mutrate.chr*
do
awk {'print $1"\t"$2-1"\t"$3"\t"$6'} $i > $i.temp
mv $i.temp $i
done