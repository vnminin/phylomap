\name{sumstatMCMCbf}
\alias{sumstatMCMCbf}
\title{Phylogenetic Stochastic Mapping for a 2x2 unknown rate matrix (without Matrix Exponentiation)}
\description{
  Creates an MCMC chain of substitution histories conditional on a phylogenetic tree, and observed states at the tips of the tree. Returns summary statisitics of each sampled substitution history.
}
\usage{
sumstatMCMCbf(z,Q,pid,Omega,N,prior)
}
\arguments{
  \item{z}{a phylogenetic tree in phytools format (a list). One of the list elements, "mapnames", must be added to the phylosim format. "mapnames" is a list of vectors. Each vector is composed of integers that are the names in phylosim's "maps" list. The "mapnames" list has a vector for each branch composed of the states visited along each branch (in order).}
  \item{Q}{a 2x2 rate matrix, the rows sum to zero}
  \item{pid}{a probability vector containing the prior probabilities that the root is in each state}
  \item{Omega}{a positive real number, the absolute value of each diagonal element of Q will be less than Omega}
  \item{N}{a positive integer: the number of substitution histories to be sampled}
  \item{prior}{a vector of length 4. The first two elements are the shape and rate parameters for the prior on Q[1,2]. The second two elements are the shape and rate parameters for the prior on Q[2,1].}
}
\details{}
\value{sumstatMCMCbf returns a matrix, each row contains summary statistics of one of the N sampled substitution histories. The first column contains the amount of time spent in state 1 over the entire tree. The second column contains the amount of time spent in state 2 over the entire tree. The third column contains the number of virtual transitions from state 1 to state 1. The fourth column contains the number of real transitions from state 1 to state 2. The fifth column contains the number of real transitions from state 2 to state 1. The sixth column contains the number of virtual transitions from state 2 to state 2. The seventh column contains the transition rate from state 1 to state 2, l01 (also called Q[1,2]). The eighth column contains the transition rate from state 2 to state 1, l10 (also called Q[2,1]). The nineth column contains the state of the root (starting with state 0, sorry for the lack of consistency). 
}
