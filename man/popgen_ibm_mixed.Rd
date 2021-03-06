% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{popgen_ibm_mixed}
\alias{popgen_ibm_mixed}
\title{popgen_ibm_mixed}
\usage{
popgen_ibm_mixed(p_net, ini_dist = NULL)
}
\arguments{
\item{p_net}{A popsnetwork object.}

\item{ini_dist}{Initial distribution of allele frequencies (optional). ini_dist has to be 
a list
containing a vector of node IDs (see \code{\link{popsnetwork}}) and
a matrix of allele frequencies. Note that *any* node pre-set in this
way will effectively be treated as a source and hide nodes that are further upstream (see
\code{\link{set_allele_freqs}}).}
}
\value{
A new popsnetwork object with allele frequencies set for each node.
}
\description{
Simulate spread of pathogens on the network using a (very) simple individual-based
model.
}
\details{
This function simulates the change of gene frequencies in a population
of pathogens as they spread through the transport network starting at the 
external sources (see \code{\link{popsnetwork}}). At each node founder
effects are assumed to change composition of the population. This change is 
simulated by directly drawing from the distribution of genotypes and unfected units,
respectively. Note that for numerical reasons this model will not work for small
absolute numbers of units (input and transfer rates). It is also strongly recommended
to only run this model on networks that have been generated with the "units" method
(see \code{\link{popsnetwork}}).
}
\examples{
# create network
el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(150, 100, 200))
ext <- data.frame(node=c("A", "B"), rate=c(300, 100), input=c(1000, 1000))
net <- popsnetwork(el, ext, spread_model="units")

# set allele frequencies (2 nodes, 3 alleles)
freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
ini_freqs <- list(as.factor(c("A", "C")), freqs)

# run simulation on initialized net
ini_net <- set_allele_freqs(net, ini_freqs)
popgen_ibm_mixed(ini_net)

# or we can initialize and run in one call
popgen_ibm_mixed(net, ini_freqs)
}
