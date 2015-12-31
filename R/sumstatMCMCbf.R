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


sumstatMCMCbf<-function(z,Q,pid,Omega,N,prior) {
 #if(dim(z$edge)[1]>9999) {
 # print("too many edges in the tree, the limit is 10,000")
 # return()
 #}

 nen=pruningwiseedgeorder(z)
 nodelist=makenodelist(z)
 root = myreorder(z)
 B=diag(dim(Q)[1])+Q/Omega

 ss<-maketreelistMCMCbf(z,Q,pid,B,Omega,nen,nodelist,root,N,prior)
 colnames(ss)<-c("time 0","time 1","n00","n01","n10","n11","l01","l10","root_state")
 return(ss[,c(1:9)])
}