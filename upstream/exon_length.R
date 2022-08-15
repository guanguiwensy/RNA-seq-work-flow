library(parallel)
library("GenomicFeatures")

#Assign annotation File and format.
gtf <- "Homo_sapiens.GRCh38.107.gtf" 
format <- "gtf"
#The GTF files are organized into a gene exon list by the GenomicFeatures package.
txdb <- makeTxDbFromGFF(gtf,format = format)
exons_gene <- exonsBy(txdb, by = "gene")

#Detect the number of threads, initialize the parallel package, and prepare for parallelization
clnum <- detectCores()-2 
cl <- makeCluster(getOption("cl.cores", clnum))

#Define a function exon_sum that can calculate the exon length of each gene in the list
exon_sum <- function(x){sum(width(reduce(x)))}
#Parallel calculation of exon length sums is achieved by parLapply. Use the function exon_sum.
exons_gene_len <- parLapply(cl,exons_gene,exon_sum)
exons_gene_len <- t(data.frame(exons_gene_len)) #List to data frame
#Save the gene length, separated by tabs.
write.table(exons_gene_len,"exons_gene_lens.txt",col.names = F, quote = F, sep = "\t")