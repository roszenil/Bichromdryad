### Assumption just one polyploidy rate for both
## size= Max number of chromosomes chosen by user
## theta= Vector of parameters lambda_i,mu_i,rho_i,p_ij
## tree.branch= length 
## p.0 initial vector (P(t)=p.0exp{Qt})
## Added two parameters to represent dynamics at states size+ that might be totally different from others.
build.sparse<-function(size,log.theta){
#Parameters
theta<-exp(log.theta)
l.0<-theta[1]
l.1<-theta[2]
m.0<-theta[3]
m.1<-theta[4]
r.0<-theta[5]
r.1<-theta[5]
prob.01<-theta[6]
prob.10<-theta[7]
e.0<-theta[8]
e.1<-theta[9]
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
	Q[i,i]<- -(m.1+l.1+r.1+prob.10)
	Q[i,(i-1)]<-m.1
	Q[i,(i+1)]<-l.1
	Q[i,2*C]<- r.1+Q[i,2*C]
	Q[i,(i-C)]<- prob.10
	}
	
Q[(2*C),(2*C)]<- -(e.1)
Q[(2*C),C]<-e.1
return(Q)
}
