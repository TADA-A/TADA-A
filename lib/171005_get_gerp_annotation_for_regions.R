# This script will take as input some genomic regions and will get all the bed intervals that have gerp score greater than some threshold.
# Usage:
# Rscript script.R [input_regions] [gerp_cutoff] [output_file]

# prefix for temporary files that will be deleted at the end of the pipeline
prefix <- system("date +%s", intern = TRUE) # prefix for temporary files that will be deleted at the end of the pipeline

args <- c("/media/yuwen/F/TADA-A/data/enhancers/Example_enhancers.bed.tempaa", "/media/yuwen/F/TADA-A/other_annotations/conservation_raw_file/hg19_gerp_score.bw",
          2, "/media/yuwen/F/TADA-A/data/enhancers/Example_enhancers.bed.tempaa.gt2.gerp.bed")
options(scipen=999)
library(data.table)
library(parallel)
args <- commandArgs(trailingOnly = TRUE)
print(args)
input_regions <- fread(args[1], header = FALSE, sep = "\t", stringsAsFactors = FALSE) 
input_regions$data_bins <- 1
data_bins <- input_regions$data_bins
input_regions <- split(input_regions, data_bins)


#funtion to expand windows to bases
window_expansion <- function(table_row){
  start <- seq(as.numeric(table_row[2]),as.numeric(table_row[3])-1)
  data.frame(table_row[1], start, start+1)
}

options(warn=-1)
cl <- makeCluster(6)

input_regions <- rbindlist(parApply(cl, input_regions[[1]], 1, window_expansion))

stopCluster(cl)
options(warn = 0)

input_regions <- data.table(input_regions, seq(1, dim(input_regions)[1]))


input_regions[[2]] <- as.character(input_regions[[2]])
input_regions[[3]] <- as.character(input_regions[[3]])
input_regions[[4]] <- as.character(input_regions[[4]])

fwrite(input_regions, paste(args[1], prefix, "temp", sep = "_"), col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
command <- paste("../../external_tools/bigWigAverageOverBed ", args[2], " ",paste(args[1], prefix, "temp", sep = "_"),
" ", paste(args[1], prefix, "temp.gerp", sep = "_"), sep = "") 
system(command)
command <- paste("awk '$4 > ", args[3], "' ", paste(args[1], prefix, "temp.gerp", sep = "_"), "| awk '{print $1}' > ", paste(args[1], prefix, "temp.gerp2", sep = "_"), sep = "")
system(command)

temp <- fread(paste(args[1], prefix, "temp.gerp2", sep = "_"), header = FALSE, sep = "\t", stringsAsFactors = FALSE)

input_regions <- input_regions[is.element(V2, temp[[1]])]

fwrite(input_regions[,1:3], args[4], col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)

command <- paste("rm ", paste(args[1], prefix, "temp*", sep = "_"), sep = "")
system(command)
quit("no")