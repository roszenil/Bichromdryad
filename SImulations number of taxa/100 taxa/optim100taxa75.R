library( "ape" )
library( "geiger" )
library( "expm" )
library( "nloptr" )
source( "masternegloglikeeps1.R" )
source( "Qmatrixwoodherb2.R" )
source("Pruning2.R")
sim.tree<-read.tree("tree100taxa75.txt") 
sim.chrom<-read.table("chrom100taxa75.txt", header=FALSE) 
last.state=50 
x.0<- log(c(0.12, 0.001, 0.25, 0.002,0.036, 0.006, 0.04,0.02, 1.792317852, 1.57e-14))
p.0<-rep(1,2*(last.state+1))/(2*(last.state+1)) 
results<-rep(0,11) 
my.options<-list("algorithm"= "NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,"maxtime"=170000000, "maxeval"=1000) 
mle<-nloptr(x0=x.0,eval_f=negloglikelihood.wh,opts=my.options,bichrom.phy=sim.tree, bichrom.data=sim.chrom,max.chromosome=last.state,pi.0=p.0) 
print(mle) 
results[1:10]<-mle$solution 
results[11]<-mle$objective 
write.table(results,file="globalmax100taxa75.csv",sep=",") 
