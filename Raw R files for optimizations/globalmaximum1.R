source( "masternegloglikeeps1.R" )
source("eudicottree.R" )
source( "Qmatrixwoodherb2.R" )
source("Pruning2.R")
library( "expm" )
library("nloptr")
bichrom.dataset<-read.table( "eudicotvals.txt",header=FALSE,sep=",",stringsAsFactors=FALSE) 
last.state=50 

x.0<- log(c(0.100228746,0.002281069,0.22588321,0.00275355,0.036806304,0.006492386,0.03854713,0.022055368,882.4422402,4.0684168))
p.0<-rep(1,2*(last.state+1))/(2*(last.state+1))
results<-rep(0,11)
my.options<-list("algorithm"="NLOPT_LN_SBPLX","ftol_rel"=1e-08,"print_level"=1,"maxtime"=170000000, "maxeval"=1000)
mle<-nloptr(x0=x.0,eval_f=negloglikelihood.wh,opts=my.options,bichrom.phy=angiosperm.tree, bichrom.data=bichrom.dataset,max.chromosome=last.state,pi.0=p.0)
print(mle) 
results[1:10]<-mle$solution 
results[11]<-mle$objective 
write.table(results,file="globalmax1.csv",sep=",") 