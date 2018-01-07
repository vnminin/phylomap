\name{simulate_2_state_tree}
\alias{simulate_2_state_tree}
\title{Tip state simulation for a 2 state model}
\description{
  Creates a tree with tip states simulated from a 2 state model.
}
\usage{
simulate_2_state_tree(101,atree,Q,pid2)
}
\arguments{
  \item{seed}{an integer}
  \item{atree}{a phylogenetic tree in the phytools format (a list). The tree contains the element, "mapnames", which must be added to the phylosim format. "mapnames" is a list of vectors. Each vector is composed of integers that are the names in phylosim's "maps" list. The "mapnames" list has a vector for each branch composed of the states visited along each branch (in order).}
  \item{Q}{a 2x2 rate matrix, the rows sum to zero}
  \item{pid}{a probability vector containing the prior probabilities that the root is in each state}
}
\details{}
\value{simulate_2_state_tree returns a tree in phylomap format with data at the tips simulated from a 2 state model.
}
