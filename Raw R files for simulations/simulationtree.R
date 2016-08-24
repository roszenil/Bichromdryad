setwd('~/Dropbox/BiChrom/LargeQmatrix/25tree')
#source("eudicottree.R")
source( "Qmatrixwoodherb4.R" )
source("simchar2.R")

library(geiger)
library(phytools)
library(expm)

N <- 50
params<-(c(0.12, 0.001, 0.25, 0.002,0.036, 0.006, 0.04,0.02, 1.792317852, 1.57E-14)) #l.0, l.1, m.0, m.1, r.0, r.1, prob.01, prob.10, e.0, e.1

Q.bichrom<-build.sparse(size=N, theta=params)
tree.taxa<-rep(0,100)

#max(nodeHeights(angiosperm.tree)) #136.901 million years 
for(i in 1:100){
tree <- sim.bdtree(b=0.26, stop="time", t=25) # Approximatedly 110 million years
write.tree(tree,file=paste("tree25time",i,".txt", sep=""))
tree.taxa[i]<-Ntip(tree)
res <- sim.char2(tree, par <- Q.bichrom, model="discrete", root=56)
write.table(res, file=paste("chrom25time",i,".txt",sep=""), row.names=TRUE,col.names=FALSE)
}
write.table(tree.taxa ,"tree.taxa25time.txt")