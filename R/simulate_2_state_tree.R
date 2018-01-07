####
#### Given a seed, a tree, a rate matrix, and a probability vector
#### containing the prior probabilities that the root is in each 
#### state, simulates tip data from a 2-state model and returns a 
#### tree with the simulated tip states.
####

simulate_2_state_tree<-function(seed,atree,Q,pid2) {

 set.seed(seed)
 sampledata2<-sample2statehistory(atree,Q,pid2) # n01 is 21  # sampledata2[7902]
 data2<-sampledata2[c(1:7901)]
 atree2<-atree

 ### create tree with simulated tip data 
 branchcount<-length(atree2$edge.length)
 tipstates<-data2[1:3951]
 atree2$states<-tipstates
 atree2$node.states<-matrix(rep(1,branchcount*2),nrow=branchcount)
 for(j in 1:length(atree2$tip.label)) {
  row<-1
  while(j!=atree2$edge[row,2]) row<-row+1
  atree2$maps[[row]]<-rep(atree2$edge.length[row]/2,2)
  names(atree2$maps[[row]])<-c(1,tipstates[j])
  cat(j/length(atree2$tip.label),"\r")
 }
 atree2$mapnames<-list()
 for(i in 1:length(atree2$maps)) atree2$mapnames[[i]]<-as.integer(names(atree2$maps[[i]]))
 for(i in 1:branchcount) if(atree2$edge[i,2]<=3951) atree2$node.states[i,2]<-tipstates[atree2$edge[i,2]]
 atree2<-makemapnames(atree2)
 return(atree2)
}

