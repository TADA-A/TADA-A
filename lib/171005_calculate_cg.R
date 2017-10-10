library(data.table)
args <- commandArgs(trailingOnly = TRUE)
print(args)
infile <- args[1] # the bed files that needs to have CG content calculted 

#calculate CG content for 50bp, 100bp 200bp and 500bp window, 
calculate_cg <-function(seq){
  letter = as.data.frame(strsplit(c(seq),split=""))
  GC_count = length(letter[letter[,1] == "C"|letter[,1] == "c" | letter[,1] == "G" | letter[,1] == "g",1])
  GC_pct = GC_count/length(letter[,1]) 
  GC_pct
}

temp <-fread(infile, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# temp1 <- temp[1:100000,]
cg_pct <- data.frame(temp[[1]],mapply(calculate_cg,temp[[2]]))
fwrite(cg_pct, paste(infile, ".cg", sep = ""), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
quit("no")