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
makenodelist<-function(x){
 rm<-NULL
 for(i in 1:(x$Nnode-1)){
  rm[i]<-reorder(x,order="pruningwise")$edge[length(x$edge[,1])-2*i,1]
 }
 return(rm)
}
myreorder<-function(x) reorder(x,order="pruningwise")$edge[dim(x$edge)[1],1]


sumstatMCMCksmt<-function(treelist,Q,pid,Omega,N,prior) {

 z<-treelist[[1]]
 edgecount<-dim(z$edge)[1]
 treecount<-length(treelist)

 nen_m<-matrix(rep(-1,edgecount*treecount),nrow=treecount)
 for(i in 1:treecount) nen_m[i,]<-pruningwiseedgeorder(treelist[[i]])
 
 tipcount<-length(z$states)
 nodelist_m<-matrix(rep(-1,(tipcount-2)*treecount),nrow=treecount)
 for(i in 1:treecount) nodelist_m[i,]<-makenodelist(treelist[[i]])

 roots<-rep(-1,treecount)
 for(i in 1:treecount) roots[i]<-myreorder(treelist[[i]])

 B=diag(dim(Q)[1])+Q/Omega
 ss<-maketreelistMCMCksmt(treelist,Q,pid,B,Omega,nen_m,nodelist_m,roots,N,prior)
 return(ss)

}