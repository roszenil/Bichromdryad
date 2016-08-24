source( "masternegloglikereduced1.R" )
source("eudicottree.R" )
source( "Qmatrixwoodherb3.R" )
source("Pruning2.R")
library( "expm" )
library("nloptr")
bichrom.dataset<-read.table( "eudicotvals.txt",header=FALSE,sep=",",stringsAsFactors=FALSE) 
last.state=50 
x.0<- log(c(0.082394414,0.010516682,0.113450749,0.025961863,0.020513114,0.012812792,0.006351607,2939.515144,349.28132))

p.0<-rep(1,2*(last.state+1))/(2*(last.state+1))
results<-rep(0,10)
my.options<-list("algorithm"="NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,"maxtime"=170000000, "maxeval"=1000)
mle<-nloptr(x0=x.0,eval_f=negloglikelihood.wh,opts=my.options,bichrom.phy=angiosperm.tree, bichrom.data=bichrom.dataset,max.chromosome=last.state,pi.0=p.0)
print(mle) 
results[1:9]<-mle$solution 
results[10]<-mle$objective 
write.table(results,file="globalreducedmax1.csv",sep=",") 