negloglikelihood.wh <- function(log.params, bichrom.phy, bichrom.data, max.chromosome, pi.0) {
	charnum=1
	Q.wh<-build.sparse(size=max.chromosome, log.theta=log.params)
	
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
