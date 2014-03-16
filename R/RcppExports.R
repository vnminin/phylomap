test <- function(x) {
    .Call('phylomap_test', PACKAGE = 'phylomap', x)
}

maketreelistMCMC <- function(x,Q,pid,B,Omega,nen,nodelist,root,N) {
    .Call('phylomap_maketreelistMCMC', PACKAGE = 'phylomap', x,Q,pid,B,Omega,nen,nodelist,root,N)
}

SPARSEmaketreelistMCMC <- function(x,Q,pid,B,Omega,nen,nodelist,root,N) {
    .Call('phylomap_SPARSEmaketreelistMCMC', PACKAGE = 'phylomap', x,Q,pid,B,Omega,nen,nodelist,root,N)
}

maketreelistEXP <- function(x,Q,pid,nen,nodelist,root,N,lefts,rights,d) {
    .Call('phylomap_maketreelistEXP', PACKAGE = 'phylomap', x,Q,pid,nen,nodelist,root,N,lefts,rights,d)
}
