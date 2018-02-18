
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


#############################################
########## For DIC Vignette #################
#############################################

findroot<-function(x) reorder(x,order="pruningwise")$edge[dim(x$edge)[1],1]

# this function counts the number of times the state increased by exactly 1
countgains<-function(statevec) return(sum(diff(statevec)==1))
countlosses<-function(statevec) return(sum(diff(statevec)==-1))

samplethebranch<-function(startstate,branchlength,Q) {
 n<-dim(Q)[1]
 statespace<-c(1:n)
 state<-startstate
 # dab: distance along branch
 dab<-0
 segmentlengths<-NULL
 segmentstates<-c(startstate)
 i<-0
 while(dab<branchlength&i<10000) {
  seg_len<-rexp(1,rate=-Q[state,state])
  dab<-dab+seg_len
  if(dab<branchlength) {
   segmentlengths<-c(segmentlengths,seg_len)
   if(n>2)  state<-(statespace[-state])[sample(c(1:(n-1)),size=1,prob=Q[state,-state])]
   if(n==2) state<-1-(state-2)
   segmentstates<-c(segmentstates,state)
  }
  if(dab>branchlength) {
   segmentlengths<-c(segmentlengths,seg_len-(dab-branchlength))
  }
  # for fear of not stopping
  i<-i+1
 }
 returnme<-list(state,segmentlengths,segmentstates)
 names(returnme)<-c("endstate","segment_lengths","segment_states")
 return(returnme)
}

sample2statehistory<-function(tree,Q,pid) {

 nodestates<-rep(0,tree$Nnode+length(tree$states))
 n<-dim(Q)[1]
 rootnode<-findroot(tree)
 #sample the root state
 rootstate<-sample(c(1:n),size=1,prob=pid)
 #record the root state
 nodestates[rootnode]=rootstate

 edgemat<-apply(reorder(tree,order="postorder")$edge,2,rev)
 edgelengths<-rev(reorder(tree,order="postorder")$edge.length)

 n01<-0 
 n10<-0
 t0<-0
 t1<-0

 for(i in 1:length(edgemat[,1])) {
  ps<-nodestates[edgemat[i,1]]
  branch<-samplethebranch(ps,edgelengths[i],Q)
   if(length(branch$segment_lengths)!=length(branch$segment_states)) {
    print("branch sampler yields inconsistent results")
    return
   }
  nodestates[edgemat[i,2]]<-branch$endstate
  n01<-n01+countgains(branch$segment_states)
  n10<-n10+countlosses(branch$segment_states)
  t0<-t0+sum(branch$segment_lengths[branch$segment_states==1])
  t1<-t1+sum(branch$segment_lengths[branch$segment_states==2])

 }

 rm<-c(nodestates,n01,n10,t0,t1)
 return(rm)
}






changemytipstates<-function(data,atree) {
 branchcount<-length(atree$edge.length)
 tipstates<-data[1:70]
 ### turn 4 states into 2 states
 tipstates<-(((data[1:70]%%2)-1)*-1)+1
 atree$states<-tipstates
 atree$node.states<-matrix(rep(1,branchcount*2),nrow=branchcount)
 for(j in 1:length(atree$tip.label)) {
  row<-1
  while(j!=atree$edge[row,2]) row<-row+1
  atree$maps[[row]]<-rep(atree$edge.length[row]/2,2)
  names(atree$maps[[row]])<-c(1,tipstates[j])
 }
 atree$mapnames<-list()
 for(i in 1:length(atree$maps)) atree$mapnames[[i]]<-as.integer(names(atree$maps[[i]]))
 for(i in 1:branchcount) if(atree$edge[i,2]<=70) atree$node.states[i,2]<-tipstates[atree$edge[i,2]]
 atree<-makemapnames(atree)
 atreeR<-atree
 return(atreeR)
}




###
### big tree DIC code
###

make2stateDICbig<-function (mat, atree, pid,ne) {
    root = myreorder(atree)
    l01 <- mean(mat[, "l01"])
    l10 <- mean(mat[, "l10"])
    Q <- matrix(c(-l01, l10, l01, -l10), nrow = 2)
    branchnumber = length(atree$maps)
    nodecount <- 2 * atree$Nnode + 1
    TP2 <- array(rep(0, 4 * branchnumber), dim = c(2, 2, branchnumber))
    for (i in 1:branchnumber) TP2[, , i] = expm(Q * atree$edge.length[i])
    n <- atree$Nnode
    edge <- atree$edge
    #ne = pruningwiseedgeorder(atree)
    PL <- matrix(rep(0, 2 * nodecount), nrow = nodecount)
    for (i in 1:length(atree$states)) {
        if (atree$states[i] == 1) 
            PL[i, ] = c(1, 0)
        if (atree$states[i] == 2) 
            PL[i, ] = c(0, 1)
    }
    S=0
    for (i in 1:n) {
     PL[edge[ne[2*i-1],1],] = (TP2[,,ne[2*i-1]] %*% PL[edge[ne[2*i-1],2],]) * (TP2[,,ne[2*i]] %*% PL[edge[ne[2*i],2],])
     S=S+log(sum(PL[edge[ne[2*i-1],1],]))
     PL[edge[ne[2*i-1],1],]=PL[edge[ne[2*i-1],1],]/sum(PL[edge[ne[2*i-1],1],])
    }
    D <- -2 * (log(as.numeric(PL[root, ] %*% pid))+S)
    pD <- mean(-2 * mat[, "log(p(y|Q))"]) - D
    DIC <- D + 2 * pD
    return(DIC)
}

make4stateDICbig<-function(mat,atree,pid,ne) {
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
  #ne=pruningwiseedgeorder(atree)

  PL<-matrix(rep(0,4*nodecount),nrow=nodecount)
  for(i in 1:length(atree$states)) {
   if(atree$states[i]==1) PL[i,]=c(1,0,1,0)
   if(atree$states[i]==2) PL[i,]=c(0,1,0,1)
  }  
  S=0
  for(i in 1:n) {
   PL[edge[ne[2*i-1],1],]=(TP2[,,ne[2*i-1]]%*%PL[edge[ne[2*i-1],2],])*(TP2[,,ne[2*i]]%*%PL[edge[ne[2*i],2],])
   S=S+log(sum(PL[edge[ne[2*i-1],1],]))
   PL[edge[ne[2*i-1],1],]=PL[edge[ne[2*i-1],1],]/sum(PL[edge[ne[2*i-1],1],])
  }

  ### D(\hat{Q})
  D<- -2*(log(as.numeric(PL[root,]%*%pid))+S)
  ### pD:
  pD<-mean(-2*mat[,"log(p(y|Q))"])-D
  ### DIC:
  DIC<-D+2*pD

  return(DIC)
}


####
#### A helper function to simplify DIC computations for the
#### squamate tree analysis. Creates summary statistics of
#### the posterior distributions for 12 different scenarios,
#### a 2 state and a 4 state analysis for three sets of 
#### priors each, once when the truth was generated by a two
#### state model and once when the truth was generated by a
#### four state model. Returns 12 DIC statistics.
####
make_12_chains<-function(seed,Qs,trees,pids,Omegas,N,priors) {

 Q<-Qs[[1]];tree<-trees[[1]];pid<-pids[[1]];Omega<-Omegas[[1]];prior<-priors[[1]]
 SIMtstr<-sumstatMCMC2sDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMtstr)<-c("t0","t1","n00","n01","n10","n11","l01","l10","root_state","log(p(y|Q))")
 saveRDS(SIMtstr,file=paste("SIM",seed,"squa2sPr.RData",sep=""))

 Q<-Qs[[1]];tree<-trees[[1]];pid<-pids[[1]];Omega<-Omegas[[1]];prior<-priors[[2]]
 SIMtstd<-sumstatMCMC2sDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMtstd)<-c("t0","t1","n00","n01","n10","n11","l01","l10","root_state","log(p(y|Q))")
 saveRDS(SIMtstd,file=paste("SIM",seed,"squa2sPd.RData",sep=""))

 Q<-Qs[[1]];tree<-trees[[1]];pid<-pids[[1]];Omega<-Omegas[[1]];prior<-priors[[3]]
 SIMtsts<-sumstatMCMC2sDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMtsts)<-c("t0","t1","n00","n01","n10","n11","l01","l10","root_state","log(p(y|Q))")
 saveRDS(SIMtsts,file=paste("SIM",seed,"squa2sPs.RData",sep=""))

 Q<-Qs[[2]];tree<-trees[[1]];pid<-pids[[2]];Omega<-Omegas[[2]];prior<-priors[[4]]
 SIMfstr<-sumstatMCMCksDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMfstr)<-c("t1","t2","t3","t4","n11","n12","n13","n14","n21","n22","n23","n24","n31","n32","n33","n34","n41","n42","n43","n44","l01","l10","k01","k10","gamma","root_state","log(p(y|Q))")
 saveRDS(SIMfstr,file=paste("SIM",seed,"squa4sPr.RData",sep=""))

 Q<-Qs[[2]];tree<-trees[[1]];pid<-pids[[2]];Omega<-Omegas[[2]];prior<-priors[[5]]
 SIMfstd<-sumstatMCMCksDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMfstd)<-c("t1","t2","t3","t4","n11","n12","n13","n14","n21","n22","n23","n24","n31","n32","n33","n34","n41","n42","n43","n44","l01","l10","k01","k10","gamma","root_state","log(p(y|Q))")
 saveRDS(SIMfstd,file=paste("SIM",seed,"squa4sPd.RData",sep=""))

 Q<-Qs[[2]];tree<-trees[[1]];pid<-pids[[2]];Omega<-Omegas[[2]];prior<-priors[[6]]
 SIMfsts<-sumstatMCMCksDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMfsts)<-c("t1","t2","t3","t4","n11","n12","n13","n14","n21","n22","n23","n24","n31","n32","n33","n34","n41","n42","n43","n44","l01","l10","k01","k10","gamma","root_state","log(p(y|Q))")
 saveRDS(SIMfsts,file=paste("SIM",seed,"squa4sPs.RData",sep=""))

 #### tree4: 112?

 Q<-Qs[[1]];tree<-trees[[2]];pid<-pids[[1]];Omega<-Omegas[[1]];prior<-priors[[1]]
 SIMtsfr<-sumstatMCMC2sDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMtsfr)<-c("t0","t1","n00","n01","n10","n11","l01","l10","root_state","log(p(y|Q))")
 saveRDS(SIMtsfr,file=paste("SIM",seed,"squa2s4Pr.RData",sep=""))

 Q<-Qs[[1]];tree<-trees[[2]];pid<-pids[[1]];Omega<-Omegas[[1]];prior<-priors[[2]]
 SIMtsfd<-sumstatMCMC2sDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMtsfd)<-c("t0","t1","n00","n01","n10","n11","l01","l10","root_state","log(p(y|Q))")
 saveRDS(SIMtsfd,file=paste("SIM",seed,"squa2s4Pd.RData",sep=""))

 Q<-Qs[[1]];tree<-trees[[2]];pid<-pids[[1]];Omega<-Omegas[[1]];prior<-priors[[3]]
 SIMtsfs<-sumstatMCMC2sDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMtsfs)<-c("t0","t1","n00","n01","n10","n11","l01","l10","root_state","log(p(y|Q))")
 saveRDS(SIMtsfs,file=paste("SIM",seed,"squa2s4Ps.RData",sep=""))

 Q<-Qs[[2]];tree<-trees[[2]];pid<-pids[[2]];Omega<-Omegas[[2]];prior<-priors[[4]]
 SIMfsfr<-sumstatMCMCksDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMfsfr)<-c("t1","t2","t3","t4","n11","n12","n13","n14","n21","n22","n23","n24","n31","n32","n33","n34","n41","n42","n43","n44","l01","l10","k01","k10","gamma","root_state","log(p(y|Q))")
 saveRDS(SIMfsfr,file=paste("SIM",seed,"squa4s4Pr.RData",sep=""))

 Q<-Qs[[2]];tree<-trees[[2]];pid<-pids[[2]];Omega<-Omegas[[2]];prior<-priors[[5]]
 SIMfsfd<-sumstatMCMCksDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMfsfd)<-c("t1","t2","t3","t4","n11","n12","n13","n14","n21","n22","n23","n24","n31","n32","n33","n34","n41","n42","n43","n44","l01","l10","k01","k10","gamma","root_state","log(p(y|Q))")
 saveRDS(SIMfsfd,file=paste("SIM",seed,"squa4s4Pd.RData",sep=""))

 Q<-Qs[[2]];tree<-trees[[2]];pid<-pids[[2]];Omega<-Omegas[[2]];prior<-priors[[6]]
 SIMfsfs<-sumstatMCMCksDICt(tree,Q,pid,Omega,N,prior)
 colnames(SIMfsfs)<-c("t1","t2","t3","t4","n11","n12","n13","n14","n21","n22","n23","n24","n31","n32","n33","n34","n41","n42","n43","n44","l01","l10","k01","k10","gamma","root_state","log(p(y|Q))")
 saveRDS(SIMfsfs,file=paste("SIM",seed,"squa4s4Ps.RData",sep=""))



 ne<-pruningwiseedgeorder(tree)

 DICtree2<-matrix(rep(0,6),nrow=2)
 DICtree4<-matrix(rep(0,6),nrow=2)

 DICtree2[1,1]<-make2stateDICbig(SIMtstr[-c(1:min(ceiling(N/10),1000)),],trees[[1]],pids[[1]],ne)
 DICtree2[2,1]<-make4stateDICbig(SIMfstr[-c(1:min(ceiling(N/10),1000)),],trees[[1]],pids[[2]],ne)

 DICtree2[1,2]<-make2stateDICbig(SIMtstd[-c(1:min(ceiling(N/10),1000)),],trees[[1]],pids[[1]],ne)
 DICtree2[2,2]<-make4stateDICbig(SIMfstd[-c(1:min(ceiling(N/10),1000)),],trees[[1]],pids[[2]],ne)

 DICtree2[1,3]<-make2stateDICbig(SIMtsts[-c(1:min(ceiling(N/10),1000)),],trees[[1]],pids[[1]],ne)
 DICtree2[2,3]<-make4stateDICbig(SIMfsts[-c(1:min(ceiling(N/10),1000)),],trees[[1]],pids[[2]],ne)


 DICtree4[1,1]<-make2stateDICbig(SIMtsfr[-c(1:min(ceiling(N/10),1000)),],trees[[2]],pids[[1]],ne)
 DICtree4[2,1]<-make4stateDICbig(SIMfsfr[-c(1:min(ceiling(N/10),1000)),],trees[[2]],pids[[2]],ne)

 DICtree4[1,2]<-make2stateDICbig(SIMtsfd[-c(1:min(ceiling(N/10),1000)),],trees[[2]],pids[[1]],ne)
 DICtree4[2,2]<-make4stateDICbig(SIMfsfd[-c(1:min(ceiling(N/10),1000)),],trees[[2]],pids[[2]],ne)

 DICtree4[1,3]<-make2stateDICbig(SIMtsfs[-c(1:min(ceiling(N/10),1000)),],trees[[2]],pids[[1]],ne)
 DICtree4[2,3]<-make4stateDICbig(SIMfsfs[-c(1:min(ceiling(N/10),1000)),],trees[[2]],pids[[2]],ne)

 saveRDS(DICtree2,file=paste("SIM",seed,"_DICtree2.RData",sep=""))
 saveRDS(DICtree4,file=paste("SIM",seed,"_DICtree4.RData",sep=""))

 return(list(DICtree2,DICtree4))

}

###
### Combines tree simulation and make_12_chains into a single step
###

simulate_estimate<-function(seed,atree,N) {
 library(expm)

 Q2<-matrix(c(-.001,.006,.001,-.006),nrow=2)
 pid2<-c(0.5,0.5)
 atree2<-simulate_2_state_tree(seed,atree,Q2,pid2)
 prior2r<-c(.55,1,.55,1)
 prior2d<-c(1,.1,1,.1)
 prior2s<-c(.01,.01,.01,.01)

 Q4<-make2sQ(.001,.001,.001,.03,16) 
 pid4<-c(0.25,0.25,0.25,0.25)
 atree4<-simulate_4_state_tree(seed,atree,Q4,pid4)
 prior4r<-c(1,10,2,10,20,2)
 prior4d<-c(1,.1,1,.1,1,.1)
 prior4s<-c(.01,.01,.01,.01,.01,.01)
 
 priors<-list(prior2r,prior2d,prior2s,prior4r,prior4d,prior4s)
 Qs<-list(Q2,Q4)
 trees<-list(atree2,atree4)
 pids<-list(pid2,pid4)
 Omegas<-list(10,10)

 DIC_list<-make_12_chains(seed,Qs,trees,pids,Omegas,N,priors)
 return(DIC_list)
}

# We prefer the model with the smallest DIC value. 
# DIC_mat is 2 by 3. First row corresponds to a 2-state DIC.
# Second row corresponds to a 4-state DIC. Columns correspond
# to restricted, diffuse, and spike priors respectively. 
correct_DIC_decision_counter<-function(DIC_mat,true_state_size) {
 if(true_state_size==2) return(as.numeric(DIC_mat[1,]<DIC_mat[2,]))
 if(true_state_size==4) return(as.numeric(DIC_mat[1,]>DIC_mat[2,]))
}

###
### Once simulate_estimate has been run 10 times with seeds 1:10,
### We combine the results into a single 2 by 3 matrix which contains
### the counts of the number of times DIC correctly choose the 2-state
### model (top row) and the number of times DIC correctly choose the
### 4 state model (bottom row). The first column corresponds to the
### restricted priors, the second column corresponds to the diffuse 
### priors, and the third column corresponds to the spike priors.
###
combine_DIC_results<-function(iterations) {
 ten_DIC_sets_tree2<-vector("list",iterations)
 for(i in 1:iterations) ten_DIC_sets_tree2[[i]]<-readRDS(paste("SIM",100+i,"_DICtree2.RData",sep=""))
 noca_tree2_list<-lapply(ten_DIC_sets_tree2,correct_DIC_decision_counter,2)
 ten_DIC_sets_tree4<-vector("list",iterations)
 for(i in 1:iterations) ten_DIC_sets_tree4[[i]]<-readRDS(paste("SIM",100+i,"_DICtree4.RData",sep=""))
 noca_tree4_list<-lapply(ten_DIC_sets_tree4,correct_DIC_decision_counter,4)
 noca_tree2_matrix <- matrix(unlist(noca_tree2_list), ncol = 3, byrow = TRUE)
 noca_tree4_matrix <- matrix(unlist(noca_tree4_list), ncol = 3, byrow = TRUE)
 correct_counts<-matrix(nrow=2,ncol=3)
 colnames(correct_counts)<-c("restricted","diffuse","spike")
 rownames(correct_counts)<-c("truth: 2 states","truth: 4 states")
 correct_counts[1,]<-apply(noca_tree2_matrix,2,sum)
 correct_counts[2,]<-apply(noca_tree4_matrix,2,sum)
 saveRDS(correct_counts,file="correct_counts.RData")
 return(correct_counts)
}


###
### corHMM AIC code
###

# We prefer the model with the smallest AIC value. 
correct_AIC_decision_counter<-function(AIC_vec,true_state_size) {
 if(true_state_size==2) return(as.numeric(AIC_vec[1]<AIC_vec[2]))
 if(true_state_size==4) return(as.numeric(AIC_vec[1]>AIC_vec[2]))
}

# 1 means corHMM's AIC result choose correctly, 0 means incorrectly
get_corHMM_AIC_result<-function(seed,atree,true_state_count) {
 library(corHMM)
 library(expm)
 if(true_state_count==2) {
  Q2 <- matrix(c(-0.001, 0.006, 0.001, -0.006), nrow = 2)
  pid2 <- c(0.5, 0.5)
  atree2 <- simulate_2_state_tree(seed, atree, Q2, pid2)
  traitmate2<-data.frame(atree2$tip.label,as.integer(atree2$state-1))
  names(traitmate2)<-c("species","T1")
  one.cat2 <- corHMM(atree2,traitmate2,rate.cat=1,node.states="marginal", ip=1)
  two.cat2 <- corHMM(atree2,traitmate2,rate.cat=2,node.states="marginal", ip=1) 
  AIC_vec<-c(one.cat2$AIC,two.cat2$AIC)
  return(correct_AIC_decision_counter(AIC_vec,2))
 }

 if(true_state_count==4) {
  Q4 <- make2sQ(0.001, 0.001, 0.001, 0.03, 16)
  pid4 <- c(0.25, 0.25, 0.25, 0.25)
  atree4 <- simulate_4_state_tree(seed, atree, Q4, pid4)
  traitmate4<-data.frame(atree4$tip.label,as.integer(atree2$state-1))
  names(traitmate4)<-c("species","T1")
  one.cat4 <- corHMM(atree4,traitmate4,rate.cat=1,node.states="marginal", ip=1)
  two.cat4 <- corHMM(atree4,traitmate4,rate.cat=2,node.states="marginal", ip=1) 
  AIC_vec<-c(one.cat4$AIC,two.cat4$AIC)
  return(correct_AIC_decision_counter(AIC_vec,4))
 }
}
