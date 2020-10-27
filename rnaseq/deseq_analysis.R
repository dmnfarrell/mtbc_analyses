# this script runs expression analytsis for rna-seq data with DESeq2
# simple case of 2 conditions
# https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html
# example tutorial - https://4va.github.io/biodatasci/r-rnaseq-airway.html

#load library
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

#set the working directory to your own folder
setwd('~/mbovis_rnaseq_jan19/')

#load phop data
#samples <- read.table("samples.csv", sep=',',header = TRUE)
#samples <-samples[1:8,]
#counts <- read.table("counts_phop.csv",header=TRUE,sep=',',row.names = 1)

#load pknh data same as above but different samples
samples <- read.table("samples.csv", sep=',',header = TRUE)
samples <- samples[9:14,]
counts <- read.table("results/counts_pknh.csv",header=TRUE,sep=',',row.names = 1)

#set the row names of the samples table
#NB: assumes they are in the same order as the counts columns
rownames(samples) <-samples$id
colnames(counts)
samples
head(counts)

#put data into a deseq object
dds <- DESeqDataSetFromMatrix(counts, samples, design = ~0+strain)
levels(dds$strain)
mycontrast=c("strain","pknH", "WT")
#run the DE
dds <- DESeq(dds)

#get results
res <- results(dds, lfcThreshold=log2(2), contrast=mycontrast )
resultsNames(dds)
mcols(res, use.names=TRUE)
summary(res)

#plots
plotMA(res)
vsd <- vst(dds) #normalised counts
write.csv(assay(vsd),'pknh_vsd.csv')

#pca
plotPCA(vsd, intgroup = "strain", ntop = 50)
#volcano plot
plot(x=res$log2FoldChange, y=-log10(res$padj), 
     xlab="log2(Fold-Change)", ylab="-log10(adjusted P-value",
     col=ifelse(res$padj<=0.05, "red", "black"), main="Volcano plot")

#get top genes
res <- res[order(res$log2FoldChange),]
res <- na.omit(res)
keep <- res[res$padj < 0.05,]
keep
tags <- as.character(sort(rownames(keep)))
id="Mb2490" 
#assay(vsd)[id,]
#plot individual counts
#plotCounts(dds, gene=id, intgroup="strain")

#heatmaps of counts for de genes
#png("heatmap_de.png")
df <- as.data.frame(colData(dds)[,c("id","strain")])

pheatmap(assay(vsd)[tags[1:20],], cluster_rows=FALSE, cluster_cols = FALSE,
         annotation_col=df["strain"], cellwidth=20, cellheight=15)

#save genes
write.csv(res,'results/de_pknh.csv')
write.csv(keep,'results/de_pknh_filtered.csv')

#code below is more specific stuff

#plot all genes in one pdf
pdf("pknh-de-plots.pdf")
for (id in tags) {
  x<-plotCounts(dds, gene=id, intgroup="strain", returnData=TRUE)  
  p <- ggplot(x, aes(x=strain, y=count, fill=strain)) +
    scale_y_log10() + 
    geom_dotplot(binaxis="y", stackdir="center") + ggtitle(id)
  print(p)
}
while (!is.null(dev.list()))  dev.off()

