\name{sumstatEXP}
\alias{sumstatEXP}
\title{Phylogenetic Stochastic Mapping with Matrix Exponentiation}
\description{
  Creates substitution histories conditional on a phylogenetic tree, a rate matrix, and observed states at the tips of the tree. Returns the sufficicent statisitics of each sampled substitution history.
}
\usage{
sumstatEXP(z,Q,pid,N)
}
\arguments{
  \item{z}{a phylogenetic tree in phytools format (a list). One of the list elements, "mapnames", must be added to the phylosim format. "mapnames" is a list of vectors. Each vector is composed of integers that are the names in phylosim's "maps" list. The "mapnames" list has a vector for each branch composed of the states visited along each branch (in order).}
  \item{Q}{a rate matrix, the rows sum to zero}
  \item{pid}{a probability vector containing the prior probabilities that the root is in each state}
  \item{N}{a positive integer: the number of substitution histories to be sampled}
}
\details{}
\value{sumstatEXP returns a matrix, each row contains the sufficient statistics of one of the N sampled substitution histories. Let there be n states total (Q and B are n by n matrices). The first n elements of each row contain the amount of time spent over the entire tree in each state (the amount of time spent in state 1 is in the first column, the amount of time spent in state 2 is in the second column, ..., the amount of time spent in state n is found in the nth column. The remaining columns contain the number of transitions between each pair of states. We start with the number of transitions from state 1 to state 2 (found in column n+1). Column n+2 contains the number of transitions from state 1 to state 3. Column n+(n-1) contains the number of transitions from state 1 to state n. Column n+(n-1)+1 contains the number of transitions from state 2 to state 1. Column n+(n-1)+2 contains the number of transitions from state 2 to state 3. Column n+(n-1)+(n-1) contains the number of transitions from state 2 to state n. Column n+(n-1)*(n-1) contains the number of transitions from state n-1 to state n. Column n+(n-1)*(n-1)+1 contains the number of transitions from state n to state 1. Column n+n*(n-1) contains the number of transitions from state n to state n-1.
}
