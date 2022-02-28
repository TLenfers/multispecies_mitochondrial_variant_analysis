#library(dplyr)
#library(stringr)
#library(gplots)
#library(ggplot2)
#library(plotrix)
#library(pheatmap)

#library(RColorBrewer)
#library(viridis)
#library(stats)


in_file <- snakemake@input[[1]]
out_file_1 <- snakemake@output[[1]]
out_file_2 <- snakemake@output[[2]]
out_file_3 <- snakemake@output[[3]]
out_file_4 <- snakemake@output[[4]]


variants <- read.table(in_file)
first_line <- readLines(in_file)
first_line <- strsplit(first_line,"\t")[[1]]
colnames(variants) <- first_line

variants$chr_pos_ref_alt <- paste(variants$`#CHROM`,variants$POS,variants$REF,variants$ALT,sep="_")
variants[,c("INFO","ID","QUAL","FILTER","FORMAT","#CHROM","POS","REF","ALT")] <- NULL


rownames(variants) <- variants$chr_pos_ref_alt
variants$chr_pos_ref_alt <- NULL
variants$N_CHR<- NULL

variants[variants==".:."]<-510

alt_variants <- variants

for(i in 1:ncol(variants)){
  variants[,i] <-gsub("1:","",variants[,i])
  variants[,i] <- gsub("\\d*,","",variants[,i])
  variants[,i] <- as.numeric(variants[,i])
}

for(i in 1:ncol(alt_variants)){
  variants[,i] <-gsub("1:","",variants[,i])
  variants[,i] <- gsub(",\\d*","",variants[,i])
  variants[,i] <- as.numeric(variants[,i])
}

names_variants <- colnames(variants)
for(i in 1:length(names_variants)){
  names_variants[i] <- gsub("\\w*\\/","",names_variants[i])
  names_variants[i] <- gsub(".bam","",names_variants[i])
}
colnames(variants) <- names_variants
colnames(alt_variants) <- names_variants



var_names <- colnames(variants)
variants <- as.matrix(variants)

my.breaks <- c(seq(0,255, by=25)) 
my.colors <- c(colorRampPalette(colors = c("red","grey","white"))(length(my.breaks)))

# heatmap with dendograms only for samples
pdf(out_file_1,         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    paper = "A4")     
pheatmap(variants,
         cluster_rows = FALSE,
         legend=TRUE,
         cellwidth = 7,
         cellheight = 1,
         fontsize_row = 1.2,
         fontsize_col = 7,
         fontsize = 7,
         color=my.colors,
         breaks=my.breaks,
         )
dev.off()

# heatmap with dendogramm for variants and dendogramm for samples

pdf(out_file_2,         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    paper = "A4")     
out <- pheatmap(variants,
         cluster_rows =TRUE,
         legend=TRUE,cellwidth = 7,
         cellheight = 1,
         fontsize_row = 1.2,
         fontsize_col = 7,
         fontsize = 7,
         color=my.colors,
         breaks=my.breaks
         )
dev.off()

############################################################################################
# heatmap for alternative variants
############################################################################################

var_names <- colnames(alt_variants)
variants <- as.matrix(alt_variants)

my.breaks <- c(seq(0,255, by=25))
my.colors <- c(colorRampPalette(colors = c("red","grey","white"))(length(my.breaks)))

# heatmap with dendograms only for samples
pdf(out_file_3,         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    paper = "A4")
pheatmap(alt_variants,
         cluster_rows = FALSE,
         legend=TRUE,
         cellwidth = 7,
         cellheight = 1,
         fontsize_row = 1.2,
         fontsize_col = 7,
         fontsize = 7,
         color=my.colors,
         breaks=my.breaks,
)
dev.off()

# heatmap with dendogramm for variants and dendogramm for samples

pdf(out_file_4,         # File name
    width = 8, height = 7, # Width and height in inches
    bg = "white",          # Background color
    paper = "A4")
out <- pheatmap(alt_variants,
                cluster_rows =TRUE,
                legend=TRUE,cellwidth = 7,
                cellheight = 1,
                fontsize_row = 1.2,
                fontsize_col = 7,
                fontsize = 7,
                color=my.colors,
                breaks=my.breaks
)
dev.off()
