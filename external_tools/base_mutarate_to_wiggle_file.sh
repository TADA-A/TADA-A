touch $1.wiggle
echo "variableStep chrom=chr1" >> $1.wiggle
awk '$1 == "chr1"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr2" >> $1.wiggle
awk '$1 == "chr2"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr3" >> $1.wiggle
awk '$1 == "chr3"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr4" >> $1.wiggle
awk '$1 == "chr4"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr5" >> $1.wiggle
awk '$1 == "chr5"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr6" >> $1.wiggle
awk '$1 == "chr6"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr7" >> $1.wiggle
awk '$1 == "chr7"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr8" >> $1.wiggle
awk '$1 == "chr8"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr9" >> $1.wiggle
awk '$1 == "chr9"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr10" >> $1.wiggle
awk '$1 == "chr10"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr11" >> $1.wiggle
awk '$1 == "chr11"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr12" >> $1.wiggle
awk '$1 == "chr12"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr13" >> $1.wiggle
awk '$1 == "chr13"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr14" >> $1.wiggle
awk '$1 == "chr14"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr15" >> $1.wiggle
awk '$1 == "chr15"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr16" >> $1.wiggle
awk '$1 == "chr16"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr17" >> $1.wiggle
awk '$1 == "chr17"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr18" >> $1.wiggle
awk '$1 == "chr18"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr19" >> $1.wiggle
awk '$1 == "chr19"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr20" >> $1.wiggle
awk '$1 == "chr20"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr21" >> $1.wiggle
awk '$1 == "chr21"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chr22" >> $1.wiggle
awk '$1 == "chr22"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chrX" >> $1.wiggle
awk '$1 == "chrX"' $1 | awk {'print $3" "$6'} >> $1.wiggle

echo "variableStep chrom=chrY" >> $1.wiggle
awk '$1 == "chrY"' $1 | awk {'print $3" "$6'} >> $1.wiggle

