% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{descendants}
\alias{descendants}
\title{descendants}
\usage{
descendants(edgelist, node)
}
\arguments{
\item{edgelist}{The graph as an edgelist (a data frame consisting of a from and 
a to column).}

\item{node}{The node to find the descendants of.}
}
\value{
An edgelist containing the complete downstream graph of the given node.
}
\description{
get all descendants of a particular node in a graph
}
\details{
This function will return the graph reachable from a given node in a
network. If no node argument is given the function is applied to each node in the
network.
}
\examples{
el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
# get descendants of specific nodes
descendants(el, "A")
descendants(el, c("A", "B"))
# apply to all nodes in the network
descendants(el)
}
