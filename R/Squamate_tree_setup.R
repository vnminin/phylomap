####
#### This code is used to convert the squamate tree into a form that
#### is compatible with phylomap. You don't need to run it as the 
#### output can be found within the extdata folder of the phylomap 
#### R package. The output is named phylomap_compatible_squamate_tree.RData.
####

make_squamate_tree<-function() {

library(phylomap)
library(ape)
squamatetree = read.tree(file=system.file("extdata/Squamate/squamate.phy",package="phylomap"))
squamateTraits = read.csv(file=system.file("extdata/Squamate/squamate_tipdata.csv",package="phylomap"),header=TRUE)

### the tree has 4162 tips
### the trait data has 8007 species
### 12 species have tip data coded as "0/1"
### 11 out of the 12 species coded as "0/1" are found at the tips of the tree
### 200 tip names are not found in the trait data csv file

### solution: drop all 211 tips as was done in Pyron 2014

checktip<-function(atree,traits,index) {
 rm<-is.element(atree$tip.label[index],traits[,1])
 return(rm)
}

### delete 200 tips not found in the trait data csv file

atree<-squamatetree
i<-1
while(i<length(atree$tip.label)) {
 if(checktip(atree,squamateTraits,i)) {
  i<-i+1
 } else {
  atree<-drop.tip(atree,i)
 }
}

### delete 11 tips with trait data labeled as "0/1"

dmvec<-NULL
for(i in 1:length(atree$tip.label)) if(squamateTraits[squamateTraits[,1]==atree$tip.label[i],2]=="0/1") dmvec<-c(dmvec,i)
for(i in rev(dmvec)) atree<-drop.tip(atree,i)

### 
### set up atree to be phylomap compatible
###

tipNames = as.character(squamateTraits[,1])
squamateTraitVec = squamateTraits[,2]
names(squamateTraitVec) = tipNames

branchcount<-length(atree$edge.length)
atree$maps<-list()
for(i in 1:length(atree$edge.length)) {
 atree$maps[[i]]<-rep(atree$edge.length[i]/100,100)
 names(atree$maps[[i]])<-c(1,1)
}
tipstates<-rep(-10,length(atree$tip.label))
for(i in 1:length(atree$tip.label)) {
 for(j in 1:length(squamateTraitVec)) {
  if(atree$tip.label[i]==names(squamateTraitVec)[j]) {
   if(squamateTraitVec[j]=="0") tipstates[i]<-1
   if(squamateTraitVec[j]=="1") tipstates[i]<-2
  }
 }
 cat(i,"\r")
}
atree$states<-tipstates
atree$node.states<-matrix(rep(1,branchcount*2),nrow=branchcount)
for(j in 1:length(atree$tip.label)) {
 row<-1
 while(j!=atree$edge[row,2]) row<-row+1
 atree$maps[[row]]<-rep(atree$edge.length[row]/100,100)
 names(atree$maps[[row]])<-c(rep(1,99),tipstates[j])
}
for(j in 1:branchcount) if(atree$edge[j,2]>length(atree$tip.label)) names(atree$maps[[j]])<-rep(1,100)
atree$mapnames<-list()
for(i in 1:length(atree$maps)) atree$mapnames[[i]]<-as.integer(names(atree$maps[[i]]))
for(i in 1:branchcount) if(atree$edge[i,2]<=length(atree$tip.label)) atree$node.states[i,2]<-tipstates[atree$edge[i,2]]
atree<-makemapnames(atree)

### Save the phylomap compatible tree
saveRDS(atree,file="phylomap_compatible_squamate_tree.RData")
#atree<-readRDS("phylomap_compatible_squamate_tree.RData")

return(atree)
}