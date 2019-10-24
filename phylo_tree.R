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
S <- region_hits[c('pro_pknh1','pro_pknh2','pknh2_sensor')]
S[S>0] = 1
states <- apply( S , 1 , paste , collapse = "-" )
states <- factor(states)

### read ANI data

anim <- read.table('ani_matrix.csv',sep=',',header=TRUE)
#rownames(anim) <- paste0(1:nrow(anim), "_", anim$query)
#anim <- anim[1:40,]
rownames(anim) <- anim$query
species <- anim['species']

# get distance matrix
dm <- dist(anim)
njtree  <- nj(dm)

# root tree
rooted <- root(njtree,"GCA_000328845")
#rooted <- root(njtree,"GCA_002865895")
rooted <- multi2di(rooted, tips=TRUE)
rooted$edge.length[rooted$edge.length <= 0] <- 1e-8

# ace function for calculating parsimony using pknh blast results
# we provide a vector containing the states - rd900 variants and 
# the phylogeny from species similarity

fit <- ace(states, rooted, type="discrete", model='ARD')

### plot tree ###

png("nj_tree_states.png",res=150, width=1500, height=1200)
#pdf("nj_tree.pdf", width=40, height=10)

# Plot the phylogeny
plot(rooted, label.offset=0.01, edge.width=1.5, show.tip.label=FALSE)
#plotTree(rooted,lwd=1,fsize=0.6,offset=1)
#title(main="Mycobacterium Species Tree")
cols<-setNames(c('red','blue','green2','orange','purple','yellow','pink'),levels(states))
tiplabels(pie=to.matrix(states[rooted$tip.label], levels(states)),
          piecol=cols, cex=0.3,col=cols)
labels <- species[rooted$tip.label,]
tiplabels(text=labels,bg=rgb(0,0,0,0),offset=.1,adj=0,frame='n')

add.simmap.legend(colors=cols,prompt=FALSE,x=0.4*par()$usr[2],
                  y=20,fsize=1.2)

#dotTree(rooted,S,standardize=TRUE,length=6,show.tip.label=FALSE)
#phylo.heatmap(rooted,S,standardize=TRUE)

#plotBranchbyTrait(rooted,states[rooted$tip.label], mode="tips", palette="heat.colors")

#nodelabels(node=1:rooted$Nnode+Ntip(rooted),
#           pie=fit$lik.anc,piecol=cols,cex=0.05)
par( mar=c(5,7,4,2), oma=c(3,1,2,2))
dev.off()
