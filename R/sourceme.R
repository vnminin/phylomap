
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