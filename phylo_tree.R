#build phylogenetic tree for MTBC genomes

library(seqinr)
library(phangorn)
library("ape")

### read data
setwd('~/gitprojects/rd900')
#anim <- read.table('ANIm_out/ANIm_percentage_identity.tab',sep='\t',header=TRUE,row.names = 1)
anim <- read.table('anim_matrix.csv',sep=',',header=TRUE)
#rownames(anim) <- paste0(1:nrow(anim), "_", anim$query)
rownames(anim) <- anim$query

### get distance matrix
dm <- dist(anim)
njtree  <- nj(dm)
plot(njtree)

# root tree
#rooted <- root(njtree,"18_H37Rv")
rooted <- root(njtree,"GCA_000195955")
rooted <- multi2di(rooted)

### plot tree ###

png("nj_tree.png",res=150, width=800, height=1200)
#tiff(file = "temp.tiff", width = 2400, height = 1400, units = "px", res = 300)

tipcol <- rep('black', length(njtree$tip.label))
# make a vector with our list of species
cats <- unique(anim$query)

# make a vector of color we want:
colorsList <-c("orange", "darkolivegreen4", "red", "orange", "turquoise", "red",
               "purple", "darkgreen", "blue","darkblue")
for(i in 1:length(cats)){
  tipcol[grep(cats[i], njtree$tip.label)] <- colorsList[i]
}

# Plotting the phylogeny
plot(rooted, label.offset=0.01, edge.width=1.5, #type='fan',
     main="Mycobacterium Species Tree", 
     tip.color=tipcol)

# Add coloured shapes for tips
#tiplabels(pch=19, 
#          col=ifelse(grepl(njtree$tip.label, pattern="mtb|CDC|H37"), "red", "blue"), 
#          cex=1)

# Set the margins
#par(mar=c(1,1,1,0))

dev.off()

#get region hits from blast analysis
region_hits = read.table('rd900_region_hits.csv',sep=',',header=TRUE)
rownames(region_hits) <- region_hits$id
# convert to discrete states
states <- region_hits[c('pknh1','pknh1.proregion','pknh2_sensor','tbd2')]
states[states>0] = 1
states <- apply( states , 1 , paste , collapse = "-" )

# ace function for calculating parsimony using pknh blast results
# we provide a vector containing the states - rd900 variants and 
# the phylogeny from species similarity  

fit <- ace(states, rooted, type="discrete")
?ace

