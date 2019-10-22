# build phylogenetic tree for MTBC genomes
# see
# http://phytools.org/mexico2018/ex/12/Plotting-methods.html
# http://blog.phytools.org/

library(seqinr)
#library(phangorn)
library("ape")
library("phytools")

setwd('~/gitprojects/rd900')

#get region hits from blast analysis

region_hits = read.table('rd900_region_hits.csv',sep=',',header=TRUE)
rownames(region_hits) <- region_hits$id
# convert to discrete states
states <- region_hits[c('pknh1','pknh1.proregion','pknh2_sensor','tbd2')]
states[states>0] = 1
states <- apply( states , 1 , paste , collapse = "-" )
states <- factor(states)

### read ANI data

#anim <- read.table('ANIm_out/ANIm_percentage_identity.tab',sep='\t',header=TRUE,row.names = 1)
anim <- read.table('ani_matrix.csv',sep=',',header=TRUE)
#rownames(anim) <- paste0(1:nrow(anim), "_", anim$query)
#anim <- anim[1:40,]
rownames(anim) <- anim$query
species <- anim['species']

### get distance matrix
dm <- dist(anim)
njtree  <- nj(dm)
plot(njtree,type='fan')

# root tree
#rooted <- root(njtree,"GCA_000195955")
rooted <- root(njtree,"GCA_002865895")
rooted <- multi2di(rooted, tips=TRUE)

### plot tree ###

png("nj_tree.png",res=150, width=800, height=1200)
#tiff(file = "temp.tiff", width = 2400, height = 1400, units = "px", res = 300)

# Plot the phylogeny
plot(rooted, label.offset=0.01, edge.width=1.5, #type='fan',
     main="Mycobacterium Species Tree", 
     tip.color=tipcol)

# Add coloured shapes for tips
#tiplabels(pch=19, 
#          col=ifelse(grepl(njtree$tip.label, pattern="mtb|CDC|H37"), "red", "blue"), 
#          cex=1)

plotTree(rooted,lwd=1,fsize=0.6,offset=1)
#title(main="Mycobacterium Species Tree")
cols<-setNames(c('red','blue','green2','orange','purple'),levels(states))
tiplabels(pie=to.matrix(states[rooted$tip.label], levels(states)),
          piecol=cols, cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.8*par()$usr[2],
                  y=60,fsize=1.2)
#plotBranchbyTrait(rooted,states[rooted$tip.label], mode="tips", palette="heat.colors")

# tip colors
tipcol <- rep('black', length(njtree$tip.label))
# make a vector with our list of species
labels <- species[rooted$tip.label,]
cats=unique(labels)
# # make a vector of color we want:
speciescolors <-c("orange", "darkolivegreen4", "red", "orange", "turquoise", "red",
                "purple", "darkgreen", "blue","darkblue")
for(i in 1:length(labels)){
   tipcol[grep(cats, labels[i])] <- speciescolors[i]
}

dev.off()

# ace function for calculating parsimony using pknh blast results
# we provide a vector containing the states - rd900 variants and 
# the phylogeny from species similarity

fit <- ace(states, rooted, type="discrete")

