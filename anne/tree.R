
setwd('/home/damien/anne')
X <- read.table('mat.csv',sep=',',header=TRUE,row.names=1)
dm <- dist(X, method="manhattan")
cl = hclust(dm)
cl
plot(cl)
hcd <- as.dendrogram(cl)

# Default plot
plot(hcd, horiz = TRUE)

library("ape")
png("tree.png",height=1200,width=1200,pointsize=20)
plot(as.phylo(cl), cex = 0.6, edge.width = 2, label.offset = 3)
dev.off()


