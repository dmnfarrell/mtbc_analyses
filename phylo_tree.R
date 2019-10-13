#build phylo tree for acinet

library(seqinr)
library(phangorn)
library("ape")

setwd('~/gitprojects/rd900')
?read.table
anim <- read.table('ANIm_out/ANIm_percentage_identity.tab',sep='\t',header=TRUE,row.names = 1)
dm <- dist(anim)
njtree  <- nj(dm)
plot(njtree)

#mymat  <- as.matrix.alignment(aln)
myboot <- boot.phylo(njtree, anim, makemytree)
rooted <- root(njtree,"GCA_001544815.1")
nodelabels(myboot,cex=0.7)          # plot the bootstrap values
njtree$node.label <- myboot   # make the bootstrap values be the node labels

plot(rooted, main = "Neighbor Joining")

###

phydata <- read.phyDat("roary_ecoli/core_gene_alignment.aln", format = "fasta")
fit = pml(rooted, phydata)
fitJC <- optim.pml(fit, TRUE)
plot(fitJC$tree)
logLik(fitJC)

bs = bootstrap.pml(fitJC, bs=10, optNni=TRUE, control=pml.control(trace = 0))
treeBS <- plotBS(midpoint(fitJC$tree), bs, p = 50, type="p")


png("ml_tree.png",res=150, width=800)
tiff(file = "temp.tiff", width = 2400, height = 1400, units = "px", res = 300)

# Set the margins
par(mar=c(1,0,0,0))

# Plotting the phylogeny
plot(treeBS, label.offset=0.0008, edge.width=2)

# Add coloured shapes for tips
tiplabels(pch=19, 
          col=ifelse(grepl(treeBS$tip.label, pattern="ecoli"), "red", "blue"), 
          cex=2)

# Add scale
add.scale.bar(lwd=3, cex=2, xpd=TRUE)
dev.off()
