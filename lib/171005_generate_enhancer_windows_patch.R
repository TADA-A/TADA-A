library(data.table)
library(parallel)
options(scipen=999)
mutrate <- fread("Example_enhancers.bed.temp.50bp.windows.bed.temp3", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
mutrate <- data.table(site_index = mutrate$V1, mutrate = mutrate$V4)
cg <- fread("Example_enhancers.bed.temp.50bp.windows.bed.fasta.cg", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(cg) <- c("site_index", "GC_content")
windows <- fread("Example_enhancers.bed.temp.50bp.windows.bed.temp", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
windows$site_index <- seq(1, length(windows[[1]]))
windows$genename <- unlist(strsplit(windows$V4, split = "_"))[seq(1, 2 * length(windows[[1]]),2)]
windows <- data.table(chr = windows$V1, start = windows$V2, end = windows$V3, site_index = windows$site_index, genename = windows$genename)
windows <- mutrate[windows, on = "site_index"]
windows <- cg[windows, on = "site_index"]
windows <-windows[,c(4,5,6,1,7,3,2)]
fwrite(windows, "Example_enhancer_windows.bed", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
quit("no")