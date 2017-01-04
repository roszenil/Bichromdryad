source("eudicottree.R") # tree file
library("expm") # exponential of the matrix
library("nloptr")
source("Pruning2.R") # pruning 
bichrom.dataset<-read.table( "eudicotvals.txt",header=FALSE,sep=",",stringsAsFactors=FALSE) #read dataset
last.state=50 # how big is our state
p.0<-rep(1,2*(last.state+1))/(2*(last.state+1))
########################################
## size= Max number of chromosomes chosen by user
## theta= Vector of parameters lambda_i,mu_i,rho_i,p_ij
## tree.branch= length 
## p.0 initial vector (P(t)=p.0exp{Qt})
build.sparse<-function(size,log.theta,log.prof){
#Parameters
theta<-exp(log.theta)
l.0<-theta[1]
l.1<-theta[2]
m.0<-theta[3]
m.1<-theta[4]
r.0<-exp(log.prof[1])
r.1<-theta[5]
prob.01<-exp(log.prof[2])
prob.10<-theta[6]
e.0<-theta[7]
e.1<-theta[8]
C<-size+1

#Building the sparse matrix
Q<-rep(0,4*C*C)
Q<- matrix(Q,ncol=2*C)
aux1<- floor(size/2)
aux2<- 2*C-1 ###############????????????

Q[1,1]<- -(l.0+r.0+ prob.01)
Q[1,2]<- l.0+r.0
Q[1,(C+1)]<-prob.01

for(i in 2:aux1){
	Q[i,i]<- -(m.0+l.0+r.0+prob.01)
	Q[i,(i-1)]<-m.0
	Q[i,(i+1)]<-l.0
	Q[i,(2*i)]<-r.0
	Q[i,(C+i)]<- prob.01
	}


for (i in (aux1+1):(C-1)){
	Q[i,i]<- -(m.0+l.0+r.0+prob.01)
	Q[i,(i-1)]<- m.0
	Q[i,(i+1)]<- l.0
	Q[i,C]<- r.0+Q[i,C]
	Q[i,(C+i)]<- prob.01
	}
	

Q[C,C]<- -(e.0)
Q[C,(2*C)]<-e.0

###################################################

Q[(C+1),(C+1)]<- -(l.1+r.1+prob.10)
Q[(C+1),(C+2)]<- l.1+r.1
Q[(C+1),1]<-prob.10

for(i in (C+2):(C+aux1)){
	Q[i,i]<- -(m.1+l.1+r.1+prob.10)
	Q[i,(i-1)]<-m.1
	Q[i,(i+1)]<-l.1
	Q[i,(2*i-C)]<-r.1
	Q[i,(i-C)]<- prob.10
	}


for(i in (C+aux1+1):aux2){
	Q[i,i]<- -(m.1+r.1+l.1+prob.10)
	Q[i,(i-1)]<-m.1
	Q[i,(i+1)]<-l.1
	Q[i,2*C]<- r.1+Q[i,2*C]
	Q[i,(i-C)]<- prob.10
	}
	
Q[(2*C),(2*C)]<- -(e.1)
Q[(2*C),C]<-e.1
return(Q)
}


negloglikeprofile <- function(log.params, bichrom.phy, bichrom.data, max.chromosome, pi.0, log.prof0) {
	charnum=1
	Q.wh<-build.sparse(size=max.chromosome, log.theta=log.params,log.prof=log.prof0)
	
	#Makes dev.rayDISC available:
	nb.tip <- length(bichrom.phy$tip.label)
	nb.node <- bichrom.phy$Nnode
	nl <- nrow(Q.wh)
	#Now we need to build the matrix of likelihoods to pass to dev.raydisc:
	liks <- matrix(0, nb.tip + nb.node, nl)
	#Now loop through the tips.
	for(i in 1:nb.tip){
		#The codon at a site for a species is not NA, then just put a 1 in the appropriate column.
		#Note: We add charnum+1, because the first column in the data is the species labels:
		if(!is.na(bichrom.data[i,charnum+1])){
			liks[i,bichrom.data[i,charnum+1]] <- 1
		}else{
			#If here, then the site has no data, so we treat it as ambiguous for all possible codons. Likely things might be more complicated, but this can modified later:
			liks[i,] <- 1
		}
	}
	#The result here is just the likelihood: 
	result <- pruning.like(p=NULL, phy=bichrom.phy, liks=liks, Q=Q.wh, rate=NULL, root.p=pi.0)
	return(result)
}
####Optimization of negloglikelihood for rho.h and q.hw

rhoqh.grid<-as.matrix(read.csv("grid33.csv",header=FALSE,sep=","))
long1<-dim(rhoqh.grid)[1]
results<-as.matrix(read.csv("resultsrhoqh33.csv",header=FALSE,sep=","))

#lambda_h=0.12618457
#lambda_w=0.001583056
# mu_h=0.249210759
#mu_w=0.002270911
#rho_h=0.036268349
#rho_w=0.006383748
#q_hw=0.040118644
#q_wh=0.022046124
#e_h=1.792317852
#e_w=1.57E-14

logmle<-log(c(0.12618457,0.001583056,0.249210759,0.002270911,0.006383748,0.022046124,1.792317852,1.57E-14))

my.options<-list("algorithm"="NLOPT_LN_SBPLX","ftol_rel"=1e-05,"print_level"=1,"maxtime"=86400,"maxeval"=500)
for(i in 18:long1){
	 a<-log(rhoqh.grid[i,])
	 mle.loop<-nloptr(x0=logmle,eval_f=negloglikeprofile, opts=my.options, bichrom.phy=angiosperm.tree,bichrom.data=bichrom.dataset, max.chromosome=last.state, pi.0=p.0,log.prof0=a)
	 print(mle.loop)
	  results[i,1:8]<-exp(mle.loop$solution)
	  results[i,9]<-mle.loop$objective  
	  log.x0<-mle.loop$solution
	  write.table(results,file="resultsrhoqh33.csv",sep=",",col.names=FALSE,row.names=FALSE)
 }


