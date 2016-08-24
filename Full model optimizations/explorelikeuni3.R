source( "masternegloglikeeps1.R" )
source("eudicottree.R" )
library( "expm" )
source( "Qmatrixwoodherb2.R" )
source("Pruning2.R")
bichrom.dataset<-read.table( "eudicotvals.txt",header=FALSE,sep=",",stringsAsFactors=FALSE) 
last.state=50 
uniform.samples<-read.csv("sample3.csv",header=FALSE) 
a<- as.numeric(t(uniform.samples)) 
p.0<-rep(1,2*(last.state+1))/(2*(last.state+1))
results<-rep(0,9)
mle<-try(optim(par=a,fn=negloglikelihood.wh, method= "Nelder-Mead", bichrom.phy=angiosperm.tree, bichrom.data=bichrom.dataset,max.chromosome=last.state,pi.0=p.0),silent=TRUE) 
print(mle) 
if(class(mle)=="try-error"){results<-rep(NA,9)}else{ 
results[1:10]<-exp(mle$par) 
results[11]<-mle$value} 
write.table(results,file="results3.csv",sep=",") 
