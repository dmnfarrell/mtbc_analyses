#build phylogenetic tree for MTBC genomes

library(seqinr)
library(phangorn)
library("ape")

### read data
setwd('~/gitprojects/rd900')
anim <- read.table('ANIm_out/ANIm_percentage_identity.tab',sep='\t',header=TRUE,row.names = 1)
anim <- read.table('anim_matrix.csv',sep=',',header=TRUE)
rownames(anim) <- paste0(1:nrow(anim), "_", anim$X)

### get distance matrix
dm <- dist(anim)
njtree  <- nj(dm)
plot(njtree)

### bootstrap nj tree
#estimate_tr <- function(m) nj(dist(m))
#njtree$tip.label <- 

#myboot <- boot.phylo(njtree, anim, estimate_tr)
# plot the bootstrap values
#nodelabels(myboot,cex=0.7) 
rooted <- root(njtree,"14_H37Rv")
# make the bootstrap values be the node labels
rooted$node.label <- myboot   
plot(rooted, main = "Neighbor Joining")

### plot tree ###

png("ml_tree.png",res=150, width=800)
#tiff(file = "temp.tiff", width = 2400, height = 1400, units = "px", res = 300)

# Set the margins
par(mar=c(1,0,0,0))

# Plotting the phylogeny
plot(rooted, label.offset=0.0008, edge.width=2)

# Add coloured shapes for tips
tiplabels(pch=19, 
          col=ifelse(grepl(njtree$tip.label, pattern="ecoli"), "red", "blue"), 
          cex=1)

# Add scale
add.scale.bar(lwd=3, cex=2, xpd=TRUE)
dev.off()

# ace function for calculating parsimony using pknh blast results
# we provide a vector containing the states - rd900 variants and the phylogeny from species similarity
ace(states, njtree)
