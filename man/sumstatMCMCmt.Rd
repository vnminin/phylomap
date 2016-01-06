\name{sumstatMCMCmt}
\alias{sumstatMCMCmt}
\title{Phylogenetic Stochastic Mapping for a 2x2 unknown rate matrix for a set of phylogenies (without Matrix Exponentiation)}
\description{
  Creates an MCMC chain of substitution histories conditional on a set of phylogenetic tree, and observed states at the tips of the tree. Returns summary statisitics of each sampled substitution history.
}
\usage{
sumstatMCMCmt(treelist,Q,pid,Omega,N,prior)
}
\arguments{
  \item{treelist}{a set of phylogenetic trees in the phytools format (a list). The class of treelist is multiPhylo. Each tree in treelist (for example, treelist[[1]]) contains the element, "mapnames", which must be added to the phylosim format. "mapnames" is a list of vectors. Each vector is composed of integers that are the names in phylosim's "maps" list. The "mapnames" list has a vector for each branch composed of the states visited along each branch (in order).}
  \item{Q}{a 2x2 rate matrix, the rows sum to zero}
  \item{pid}{a probability vector containing the prior probabilities that the root is in each state}
  \item{Omega}{a positive real number, the absolute value of each diagonal element of Q will be less than Omega}
  \item{N}{a positive integer: the number of substitution histories to be sampled}
  \item{prior}{a vector of length 4. The first two elements are the shape and rate parameters for the prior on Q[1,2]. The second two elements are the shape and rate parameters for the prior on Q[2,1].}
}
\details{}
\value{sumstatMCMCmt returns a matrix, each row contains summary statistics of one of the N sampled substitution histories. The first column contains the amount of time spent in state 0 over the entire tree. The second column contains the amount of time spent in state 1 over the entire tree. The third column contains the number of virtual transitions from state 0 to state 0. The fourth column contains the number of real transitions from state 0 to state 1. The fifth column contains the number of real transitions from state 1 to state 0. The sixth column contains the number of virtual transitions from state 1 to state 1. The seventh column contains the transition rate from state 0 to state 1, l01 (also called Q[1,2]). The eighth column contains the transition rate from state 1 to state 0, l10 (also called Q[2,1]). The nineth column contains an integer indicating which tree was used to sample each row's substitution history (starting with tree 0 in the tree list). 
}
