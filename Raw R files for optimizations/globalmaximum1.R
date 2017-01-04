source( "masternegloglikeeps1.R" )
source("eudicottree.R" )
source( "Qmatrixwoodherb2.R" )
source("Pruning2.R")
library( "expm" )
library("nloptr")
bichrom.dataset<-read.table( "eudicotvals.txt",header=FALSE,sep=",",stringsAsFactors=FALSE) 
last.state=50 

x.0<-log(c(0.125975741,0.002013321,	0.263220692,0.001995657,0.034083811,0.008359046,0.041999762,0.02156642,40.68839451, 3.925946466))
p.0<-rep(1,2*(last.state+1))/(2*(last.state+1))
results<-rep(0,11)
my.options<-list("algorithm"="NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,"maxtime"=170000000, "maxeval"=1000)
mle<-nloptr(x0=x.0,eval_f=negloglikelihood.wh,opts=my.options,bichrom.phy=angiosperm.tree, bichrom.data=bichrom.dataset,max.chromosome=last.state,pi.0=p.0)
print(mle) 
results[1:10]<-mle$solution 
results[11]<-mle$objective 
write.table(results,file="globalmax1.csv",sep=",") 