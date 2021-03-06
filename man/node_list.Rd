% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{node_list}
\alias{node_list}
\title{node_list}
\usage{
node_list(p_net, as_string = FALSE)
}
\arguments{
\item{p_net}{A PopsNet object.}

\item{as_string}{Whether to return the list of nodes as string vector. If this is FALSE
(the default) a factor or plain integer vector (depending on how the net was constructed)
will be returned.}
}
\value{
A dataframe with id and rate_infected.
}
\description{
Get a list of nodes in a dataframe.
}
\details{
This function returns a list of the nodes in the network as well as the amount of
infected material they contain.
}
\examples{
# create network
el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
net <- popsnetwork(el, ext)
node_list(net)
}
