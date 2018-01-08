####
#### Given a seed, a tree, a rate matrix, and a probability vector
#### containing the prior probabilities that the root is in each 
#### state, simulates tip data from a 4-state model and returns a 
#### tree with the simulated tip states mapped to a 2-state model.
####
simulate_4_state_tree<-function(seed,atree,Q,pid4) {

 set.seed(seed)
 sampledata4<-sample2statehistory(atree,Q,pid4) # n01 is 112   # sampledata4[7902]
 data4<-sampledata4[c(1:7901)]

 atree4<-atree

 ### create tree with simulated tip data 
 branchcount<-length(atree4$edge.length)
 tipstates<-data4[1:3951]
 ### turn 4 states into 2 states
 tipstates<-(((data4[1:3951]%%2)-1)*-1)+1
 atree4$states<-tipstates
 atree4$node.states<-matrix(rep(1,branchcount*2),nrow=branchcount)
 for(j in 1:length(atree4$tip.label)) {
  row<-1
  while(j!=atree4$edge[row,2]) row<-row+1
  atree4$maps[[row]]<-rep(atree4$edge.length[row]/2,2)
  names(atree4$maps[[row]])<-c(1,tipstates[j])
  cat(j,"\r")
 }
 atree4$mapnames<-list()
 for(i in 1:length(atree4$maps)) atree4$mapnames[[i]]<-as.integer(names(atree4$maps[[i]]))
 for(i in 1:branchcount) if(atree4$edge[i,2]<=3951) atree4$node.states[i,2]<-tipstates[atree4$edge[i,2]]
 atree4<-makemapnames(atree4)
 return(atree4)
}

