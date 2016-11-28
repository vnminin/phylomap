
########################### going from diversitree format to phytools format: ####################################
##################################################################################################################

divtophymap<-function(phy,y){
 tool<-phy
 tool$orig<-NULL
 tool$node.label<-NULL
 tool$maps<-vector("list",dim(phy$edge)[1])
 for(i in 1:(dim(phy$edge)[1])) tool$maps[[i]]<-onemapbranch(phy,y,i)
 return(tool)
}

onemapbranch<-function(phy,y,i){
 hist<-y$history[[i]]
 hist2<-rbind(hist,c(phy$edge.length[i],hist[dim(hist)[1],2]))
 returnme<-diff(hist2[,1])
 names(returnme)<-hist[,2]
 return(returnme)
}

nodestatesmake<-function(tool){
 node.states<-tool$edge
 for(i in 1:(dim(tool$edge)[1])) node.states[i,]<-c(as.integer(names(tool$maps[[i]])[1]),as.integer(names(tool$maps[[i]])[length(tool$maps[[i]])]))
 tool$node.states<-node.states
 return(tool)
}

statesmake<-function(tool){
states<-rep(0,length(tool$tip.label))
for(i in 1:length(tool$tip.label)) states[i]<-tool$node.states[,2][tool$edge[,2]==i]
tool$states<-states
return(tool)
}

makemappededge<-function(tree,pid){
 numberofstates<-length(pid)
 returnme<-matrix(rep(0,numberofstates*length(tree$maps)),nrow=length(tree$maps))
 for(i in 1:length(tree$maps)){
  for(j in 1:numberofstates){
   returnme[i,j]<-sum(tree$maps[[i]][names(tree$maps[[i]])==j])
  }
 }
 return(returnme)
}

divtophy<-function(phy,y){
 tool<-divtophymap(phy,y)
 tool<-nodestatesmake(tool)
 tool<-statesmake(tool)
 tool$mapped.edge<-makemappededge(tool,y$states)
 return(tool)
}

makemapnames<-function(z) {
 names<-vector(mode="list",length(z$maps))
 for(i in 1:length(z$maps)) names[[i]]<-as.integer(names(z$maps[[i]]))
 z$mapnames<-names
 return(z)
}
##################################################################################################
##################################################################################################

### plotting

makecounts<-function(vec,bins) {
 rm<-NULL
 for(i in 2:(length(bins))) rm<-c(rm,length(vec[bins[i-1]<vec&vec<bins[i]]))
 return(rm)
}
makebpm<-function(mycl,bins) {
 n=length(mycl)
 names<-rep(0,length(bins)-1)
 for(i in 1:(length(bins)-1)) names[i]<-mean(c(bins[i],bins[i+1]))
 rownames<-rep("row",n)
 rm<-matrix(rep(0,n*(length(bins)-1)),nrow=n,dimnames=list(rownames,names))
 for(i in 1:n) rm[i,]<-makecounts(mycl[[i]],bins)
 names<-NULL
 return(rm)
}
multhist<-function(numjumpsEXP,numjumpsMCMC,numjumpsSPARSE) {
 l<-min(numjumpsEXP,numjumpsMCMC,numjumpsSPARSE)
 u<-max(numjumpsEXP,numjumpsMCMC,numjumpsSPARSE)
 breaks<-seq(from=floor(l)-.5,to=ceiling(u)+.5,by=1)
 barplot(makebpm(list(numjumpsEXP,numjumpsMCMC,numjumpsSPARSE),breaks),beside=TRUE,col=c("purple","black","red"),xlab="number of jumps")
 legend("topright",legend=c("EXP","MCMC","SPARSE"),col=c("purple","black","red"),lty=1,lwd=30,cex=2,bty="n")
}
##################################################################################################
##################################################################################################

### takes an ape phylogeny with state history ($maps and $mapnames) and 
### returns a phylogeny with transitions only at the midpoint of branches leading to tip nodes
crazystart1<-function(z) {
 myz<-z
 myz$node.states<-matrix(rep(1,dim(myz$node.states)[1]*dim(myz$node.states)[2]),nrow=dim(myz$node.states)[1])
 for(i in 1:dim(myz$node.states)[1]) if(myz$edge[i,2]<=length(myz$states)) myz$node.states[i,2]<-myz$states[myz$edge[i,2]]
 for(i in 1:length(myz$maps)) {
  p<-runif(1)
  myz$maps[[i]]<-c(sum(myz$maps[[i]])*p,sum(myz$maps[[i]])*(1-p))
  names(myz$maps[[i]])<-myz$node.states[i,]
 }
 for(i in 1:length(myz$mapnames)) myz$mapnames[[i]]<-as.numeric(myz$node.states[i,])
 for(i in 1:dim(myz$mapped.edge)[1]) {
  myz$mapped.edge[i,1]<-sum((myz$maps[[i]])[myz$mapnames[[i]]==1])
  myz$mapped.edge[i,2]<-sum((myz$maps[[i]])[myz$mapnames[[i]]==2])
 }
 return(myz)
}

##################################################################################################
##################################################################################################
#################                                         ########################################
################# for calculating observed data liklihood ########################################
#################                                         ########################################
##################################################################################################
##################################################################################################

### returns the node number corresponding to the root
myreorder<-function(x) reorder(x,order="pruningwise")$edge[dim(x$edge)[1],1]

### returns a vector of edge numbers: starting at the tips working up the tree visit the edges in the order specified
pruningwiseedgeorder<-function(x){
 returnme<-seq(1:dim(x$edge)[1])
 y<-reorder(x,order="pruningwise")
 for(i in 1:dim(x$edge)[1]){
  for(j in 1:dim(x$edge)[1]){
   if(x$edge[i,1]==y$edge[j,1]&x$edge[i,2]==y$edge[j,2]) returnme[j]<-i
  }
 }
 return(returnme)
}

###
### this code requires the expm library for exponentiating matrices
###

###
### 2 state 
###

make2stateDIC<-function(mat,atree,pid) {
 root = myreorder(atree)
 l01<-mean(mat[,"l01"])
 l10<-mean(mat[,"l10"])
 Q<-matrix(c(-l01,l10,l01,-l10),nrow=2)

 ### make transition probability matrices
 branchnumber=length(atree$maps) 
 nodecount<-2*atree$Nnode+1
 TP2<-array(rep(0,4*branchnumber),dim=c(2,2,branchnumber))
 for(i in 1:branchnumber) TP2[,,i] = expm(Q*atree$edge.length[i])

 ### make Partial Likelihood matrix
 n<-atree$Nnode
 edge<-atree$edge
 ne=pruningwiseedgeorder(atree)

 PL<-matrix(rep(0,2*nodecount),nrow=nodecount)
 for(i in 1:length(atree$states)) {
  if(atree$states[i]==1) PL[i,]=c(1,0)
  if(atree$states[i]==2) PL[i,]=c(0,1)
 }  

 for(i in 1:n) PL[edge[ne[2*i-1],1],]=(TP2[,,ne[2*i-1]]%*%PL[edge[ne[2*i-1],2],])*(TP2[,,ne[2*i]]%*%PL[edge[ne[2*i],2],])

 ### D(\hat{Q})
 D<- -2*log(as.numeric(PL[root,]%*%pid))

 ### pD:
 pD<-mean(-2*mat[,"log(p(y|Q))"])-D

 ### DIC:
 DIC<-D+2*pD

 return(DIC)

}

### a version of DIC that conditions on m, the number of jumps on each branch

make2stateDICm<-function(mat,atree,pid,Omega) {
 root = myreorder(atree)
 l01<-mean(mat[,"l01"])
 l10<-mean(mat[,"l10"])
 Q<-matrix(c(-l01,l10,l01,-l10),nrow=2)
 B<-diag(dim(Q)[1])+Q/Omega

 ### make transition probability matrices
 branchnumber=length(atree$maps) 
 nodecount<-2*atree$Nnode+1
 TP2<-array(rep(0,4*branchnumber),dim=c(2,2,branchnumber))

 ### median branch transition counts
 btc<-as.numeric(round(apply(mat[,-c(1:11)],2,median)))

 for(i in 1:branchnumber) TP2[,,i] = B%^%btc[i]
 #for(i in 1:branchnumber) TP2[,,i] = expm(Q*atree$edge.length[i])

 ### make Partial Likelihood matrix
 n<-atree$Nnode
 edge<-atree$edge
 ne=pruningwiseedgeorder(atree)

 PL<-matrix(rep(0,2*nodecount),nrow=nodecount)
 for(i in 1:length(atree$states)) {
  if(atree$states[i]==1) PL[i,]=c(1,0)
  if(atree$states[i]==2) PL[i,]=c(0,1)
 }  

 for(i in 1:n) PL[edge[ne[2*i-1],1],]=(TP2[,,ne[2*i-1]]%*%PL[edge[ne[2*i-1],2],])*(TP2[,,ne[2*i]]%*%PL[edge[ne[2*i],2],])

 ### D(\hat{Q})
 D<- -2*log(as.numeric(PL[root,]%*%pid))

 ### pD:
 pD<-mean(-2*mat[,"log(p(y|Q,m))"])-D

 ### DIC:
 DIC<-D+2*pD

 return(DIC)

}

###
### 4 state 
###

make2sQ<-function(l01,l10,rkappas,lkappas,gammas) {

 k<-length(rkappas)
 n<-2*k+2
 Q<-matrix(rep(0,n^2),nrow=n)
 Q[1,2]<-l01
 Q[2,1]<-l10
 for(i in 1:k) Q[2*i-1,2*i+1]<-rkappas[i]
 for(i in 1:k) Q[2*i,2*i+2]<-rkappas[i]
 for(i in 1:k) Q[2*i+1,2*i-1]<-lkappas[i]
 for(i in 1:k) Q[2*i+2,2*i]<-lkappas[i]
 for(i in 1:k) Q[2*i+1,2*i+2]<-gammas[i]*l01
 for(i in 1:k) Q[2*i+2,2*i+1]<-gammas[i]*l10
 diag(Q)<-0
 diag(Q)<-apply(-Q,1,sum)

 return(Q)
}

make4stateDIC<-function(mat,atree,pid) {
  root = myreorder(atree)
  l01<-mean(mat[,"l01"])
  l10<-mean(mat[,"l10"])
  k01<-mean(mat[,"k01"])
  k10<-mean(mat[,"k10"])
  gamma<-mean(mat[,"gamma"])
  Q<-make2sQ(l01,l10,k01,k10,gamma) 

  ### make transition probability matrices
  branchnumber=length(atree$maps) 
  nodecount<-2*atree$Nnode+1
  TP2<-array(rep(0,16*branchnumber),dim=c(4,4,branchnumber))
  for(i in 1:branchnumber) TP2[,,i] = expm(Q*atree$edge.length[i])

  ### make Partial Likelihood matrix
  n<-atree$Nnode
  edge<-atree$edge
  ne=pruningwiseedgeorder(atree)

  PL<-matrix(rep(0,4*nodecount),nrow=nodecount)
  for(i in 1:length(atree$states)) {
   if(atree$states[i]==1) PL[i,]=c(1,0,1,0)
   if(atree$states[i]==2) PL[i,]=c(0,1,0,1)
  }  

  for(i in 1:n) PL[edge[ne[2*i-1],1],]=(TP2[,,ne[2*i-1]]%*%PL[edge[ne[2*i-1],2],])*(TP2[,,ne[2*i]]%*%PL[edge[ne[2*i],2],])

  ### D(\hat{Q})
  D<- -2*log(as.numeric(PL[root,]%*%pid))
  ### pD:
  pD<-mean(-2*mat[,"log(p(y|Q))"])-D
  ### DIC:
  DIC<-D+2*pD

  return(DIC)
}


### a version of DIC that conditions on m, the number of jumps on each branch

make4stateDICm<-function(mat,atree,pid,Omega) {
 root = myreorder(atree)
 l01<-mean(mat[,"l01"])
 l10<-mean(mat[,"l10"])
 k01<-mean(mat[,"k01"])
 k10<-mean(mat[,"k10"])
 gamma<-mean(mat[,"gamma"])
 Q<-make2sQ(l01,l10,k01,k10,gamma) 
 B<-diag(dim(Q)[1])+Q/Omega

 ### make transition probability matrices
 branchnumber=length(atree$maps) 
 nodecount<-2*atree$Nnode+1
 TP2<-array(rep(0,16*branchnumber),dim=c(4,4,branchnumber))

 ### median branch transition counts
 btc<-as.numeric(round(apply(mat[,-c(1:28)],2,median)))

 for(i in 1:branchnumber) TP2[,,i] = B%^%btc[i]

 ### make Partial Likelihood matrix
 n<-atree$Nnode
 edge<-atree$edge
 ne=pruningwiseedgeorder(atree)

 PL<-matrix(rep(0,4*nodecount),nrow=nodecount)
 for(i in 1:length(atree$states)) {
  if(atree$states[i]==1) PL[i,]=c(1,0,1,0)
  if(atree$states[i]==2) PL[i,]=c(0,1,0,1)
 }  

 for(i in 1:n) PL[edge[ne[2*i-1],1],]=(TP2[,,ne[2*i-1]]%*%PL[edge[ne[2*i-1],2],])*(TP2[,,ne[2*i]]%*%PL[edge[ne[2*i],2],])

 ### D(\hat{Q})
 D<- -2*log(as.numeric(PL[root,]%*%pid))

 ### pD:
 pD<-mean(-2*mat[,"log(p(y|Q,m))"])-D

 ### DIC:
 DIC<-D+2*pD

 return(DIC)

}
