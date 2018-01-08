\name{simulate_4_state_tree}
\alias{simulate_4_state_tree}
\title{Tip state simulation for a 4 state model}
\description{
  Creates a tree with tip states simulated from a 4 state model that are mapped back to 2 states.
}
\usage{
simulate_4_state_tree(101,atree,Q,pid4)
}
\arguments{
  \item{seed}{an integer}
  \item{atree}{a phylogenetic tree in the phytools format (a list). The tree contains the element, "mapnames", which must be added to the phylosim format. "mapnames" is a list of vectors. Each vector is composed of integers that are the names in phylosim's "maps" list. The "mapnames" list has a vector for each branch composed of the states visited along each branch (in order).}
  \item{Q}{a 4x4 rate matrix, the rows sum to zero}
  \item{pid}{a probability vector containing the prior probabilities that the root is in each state}
}
\details{}
\value{simulate_4_state_tree returns a tree in phylomap format with data at the tips simulated from a 4 state model mapped back to 2 states. The first two states are mapped to state 1 and the second two states are mapped to state 2.
}
