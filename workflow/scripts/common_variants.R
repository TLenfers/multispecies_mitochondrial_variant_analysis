# read in_file lines
x <- readLines(snakemake@input[['in_file']])
chr <- 0
for (i in 1:length(x)){
  if (startsWith(x[i],"#CHROM")){
    chr <- i
  }
}
x <- x[chr:length(x)]
fileConn <- file(snakemake@output[['out']], "w")
writeLines(x,fileConn)
close(fileConn)
