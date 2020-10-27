#gene categories for pknh
library(dplyr)
library(plyr)
library(pheatmap)
library(RColorBrewer)

setwd('~/mbovis_rnaseq_jan19/')

#annotations
ann <- read.table('mbovis_categories.csv', sep=',', header = TRUE, quote='"')
ann <- na.omit(ann,cols='Mb_tag')
ann <- ann %>% distinct(ann$Mb_tag, .keep_all = TRUE)
rownames(ann) <- ann$Mb_tag

#read in DE genes table
res <- read.table("de_pknh_all.csv", sep=',', header=TRUE)
res$tag<-res$X
res <- merge(res,ann[,1:5],by.x='tag',by.y="Mb_tag",all.y=FALSE)
res$change <- cut(res$log2FoldChange, breaks = c(-5,0,5), labels=c('down','up'))
res <- res[order(res$log2FoldChange),]
write.csv(res,'pknh_genes_all.csv')

#filtered
keep <- read.table("de_pknh_filtered.csv", sep=',', header=TRUE)
keep$tag<-keep$X
head(keep)

#merge gene lists with annotation table

keep <- merge(keep,ann[,1:5],by.x='tag',by.y="Mb_tag",all.y=FALSE)
keep$change <- cut(keep$log2FoldChange, breaks = c(-5,0,5), labels=c('down','up'))
keep <- keep[order(keep$log2FoldChange),]
write.csv(keep,'pknh_genes_filtered.csv')

summ <- plyr::count(keep, c('change','category'))
summ
#qplot( data=summ, x = category, y = freq, geom="bar") + facet_wrap("change")
keep
#top de genes
f<-filter(keep,abs(log2FoldChange)>2.5)
par(mar=c(4, 15 ,5 ,12))
tags <- as.character(f$tag)

colors <- colorRampPalette( rev(brewer.pal(3, "RdBu")) )(2)
barplot(f$log2FoldChange, beside=FALSE, names.arg=f$gene, horiz=TRUE, 
        las=2, col=ifelse(f$log2FoldChange>=0,"#C11B17","#0041C2"), 
        space=.3, xlab="log2FC")

#by category
cats=levels(keep$category)[2:10]
cat="cell wall and cell processes"
f<-filter(keep, category==cat)
tags <- as.character(f$tag)
par(mar=c(4, 5 ,5 ,2))
barplot(f$log2FoldChange, beside=FALSE, names.arg=f$gene, horiz=TRUE, las=2,
        space=0, main=cat, xlab="log2FC")


#plot by category
df <- as.data.frame(colData(dds)[,c("id","strain")])
colors <- colorRampPalette( rev(brewer.pal(9, "RdBu")) )(9)

tags <- c("Mb1298C","Mb1317","Mb3646C","Mb3645C","Mb3644C","Mb1316","Mb0178","Mb0177","Mb0176","Mb0175",
          "Mb1160C","Mb1161","Mb1162","Mb2898","Mb2899","Mb2900")
genes <- c("Mb1298C","Mb1317","espA","espC","espD","Mb1316","Mb0178","Mb0177","Mb0176","Mb0175",
           "Mb1160C","prpD","prpC","mpb83","dipZ","mpb70")
x<-assay(vsd)[tags,]
rownames(x) <- genes
x

png("heatmap.png",width = 6, height = 5, units = 'in', res = 300)
pheatmap(x, cluster_rows=F, cluster_cols = FALSE, color=colors, 
         show_colnames = F, fontsize_row = 8, border_color="black",
         annotation_col=df["strain"], cellwidth=15, cellheight=8)
dev.off()

pdf("pknh-by-category.pdf")
for (cat in cats) {
  f<-filter(keep, category==cat)
  par(mar=c(4, 5 ,5 ,2))
  #p1 <- barplot(height=f$log2FoldChange, horiz=TRUE, col='lightblue', las=2, beside = TRUE,
  #        space=0, main=cat, names.arg=f$gene, xlab="log2FC")
  #print(p1)
  tags <- as.character(f$tag)
  x<-assay(vsd)[tags,]
  rownames(x) <- f$gene
  p2 <- pheatmap(x, cluster_rows=FALSE,  cluster_cols = FALSE, color=colors,
                 show_colnames = F, fontsize_row = 8, border_color="black",
                 annotation_col=df["strain"], cellwidth=20, cellheight=8, main=cat)
  print(p2)
}
while (!is.null(dev.list()))  dev.off()

