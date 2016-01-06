\name{sumstatMCMCksmt}
\alias{sumstatMCMCksmt}
\title{Phylogenetic Stochastic Mapping for a 2 state trait with an unknown rate matrix and a set of phylogenies (without Matrix Exponentiation)}
\description{
  Creates an MCMC chain of substitution histories conditional on a set of phylogenetic tree, and observed states at the tips of the tree. Returns summary statisitics of each sampled substitution history. This function uses a hidden rates model of evolution allowing for different transition rates across the trees.
}
\usage{
sumstatMCMCksmt(treelist,Q,pid,Omega,N,prior)
}
\arguments{
  \item{treelist}{a set of phylogenetic trees in the phytools format (a list). The class of treelist is multiPhylo. Each tree in treelist (for example, treelist[[1]]) contains the element, "mapnames", which must be added to the phylosim format. "mapnames" is a list of vectors. Each vector is composed of integers that are the names in phylosim's "maps" list. The "mapnames" list has a vector for each branch composed of the states visited along each branch (in order).}
  \item{Q}{an nxn rate matrix where n is an even integer strictly greater than 2, the rows sum to zero}
  \item{pid}{a probability vector containing the prior probabilities that the root is in each state}
  \item{Omega}{a positive real number, the absolute value of each diagonal element of Q will be less than Omega}
  \item{N}{a positive integer: the number of substitution histories to be sampled}
  \item{prior}{a vector of length 8. The first two elements are the shape and rate parameters for the prior on Q[1,2] (l01). The second two elements are the shape and rate parameters for the prior on Q[2,1] (l10). The third two elements are the shape and rate parameters for all the priors on the rates of transitioning between regimes (kappas). The final two elements are the shape and rate parameters for all the priors governing the multiplicative factors that distinguish the different regime transition rates (gammas).}
}
\details{}
\value{sumstatMCMCksmt returns a matrix, each row contains a summary statistic of one of the N sampled substitution histories. The first n elements of each row contain the amount of time spent over the entire tree in each state. The amount of time spent in state 1 is in the first column, the amount of time spent in state 2 is in the second column, ... , and the amount of time spent in state n is found in the nth column. The next n^2 columns contain the number of transitions between each pair of states (including virtual transitions). We start with the number of transitions from state 1 to state 1, then from state 1 to state 2, on to the number of transitions from state 1 to state n. We continue with the number of transitions from state 2 to state 1, then from state 2 to state 2, and on to the number of transitions from state 2 to state n. After n+n^2 columns we move on to reporting rate matrix parameters. In column n+n^2+1 we report the forward transition rate, Q[1,2] (l01). In column n+n^2+2 we report the backward transition rate, Q[2,1] (l10). The following (n/2)-1 columns contain regime forward transitions rates (kappa_{i(i+1)}). The following (n/2)-1 columns contain regime backward transitions rates (kappa_{(i+1)i}). The following (n/2)-1 columns contain the multiplicative factors that distinguish the different regime transition rates (gamma_i). The final column contains an integer indicating which tree was used to sample each row's substitution history (starting with tree 0 in the tree list).
}
