#ifndef DIR_NETWORK_H
#define DIR_NETWORK_H

#include "rpathsonpaths_types.h"

#include <Rcpp.h>

using namespace Rcpp;

//' @title sources
//'
//' @description Get a list of source nodes from a network.
//'
//' @details Extract a list of source nodes, i.e. nodes that do not have any
//' inputs, from a network.
//'
//' @param edge_list A dataframe containing a list of edges (see \code{\link{popsnetwork}}
//' for a description of possible formats).
//'
//' @return A list of ids of all source nodes in the network.
//'
//' @examples
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
//' sources(el)
// [[Rcpp::export]]
IntegerVector sources(const DataFrame & edge_list);


//' @title sinks 
//'
//' @description Get a list of sink (i.e. leave) nodes from a network.
//'
//' @details Extract a list of sink nodes, i.e. nodes that do not have any
//' outputs, from a network.
//'
//' @param edge_list A dataframe containing a list of edges (see \code{\link{popsnetwork}}
//' for a description of possible formats).
//'
//' @return A list of ids of all sink nodes in the network.
//'
//' @examples
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
//' sinks(el)
// [[Rcpp::export]]
IntegerVector sinks(const DataFrame & edge_list);


//' @title colour_network
//'
//' @description Identify separate sub-networks.
//'
//' @details This function identifies completely separate sub-networks in a network
//' described as an edge list.
//'
//' @param edge_list A dataframe containing a list of edges (see \code{\link{popsnetwork}}
//' for a description of possible formats).
//'
//' @return An integer vector with the sub-network id of each edge. Note that id's start at
//' 1 and are not guaranteed to be contiguous.
// [[Rcpp::export]]
IntegerVector colour_network(const DataFrame & edge_list);


//' @title cycles
//' 
//' @description Detect cycles in a network.
//' 
//' @details This function detects circular connections in a network.
//'
//' @param edge_list A dataframe containing a list of edges (see \code{\link{popsnetwork}}
//' for a description of possible formats).
//' @param record Whether to return a list of cycles.
//' 
//' @return If record is FALSE: TRUE if a cycle was found, FALSE otherwise. If record is TRUE:
//' a list of cycles (as vectors of node ids, see \code{\link{popsnetwork}}) is returned.
//'
//' @examples
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
//' cycles(el)
// [[Rcpp::export]]
SEXP cycles(const DataFrame & edge_list, bool record=false);


//' @title popsnetwork 
//'
//' @description Create a popsnetwork object.
//'
//' @details A popsnetwork object stores the nodes and edges making up a food transport
//' network and associated data describing spread of infected material on the
//' network. A network is created from a tabular description of its edges.
//'
//' @section Edge lists and node ids:
//'
//' Many functions in rpathsonpaths take edge lists as an argument. An edge list is a 
//' data frame with at least two columns. These columns can be either integer vectors
//' which will be interpreted as \emph{0-based node indices} or factors. In general where
//' applicable functions will return the same format they received. 
//'
//' Note that formats 
//' can not be mixed easily. Function arguments (and data frame columns) have to have the 
//' same type. Furthermore nodes in a \code{popsnetwork} object generated with integer
//' indices can not be referenced by name.
//'
//' In general integer indices are considerably faster than factors.
//'
//' @section Using popsnetwork:
//'
//' Generating simulated genetic data from a popsnetwork object happens in several steps:
//' \enumerate{
//' \item Create the basic network using the constructor \code{\link{popsnetwork}}.
//' \item Set initial allele frequencies for a number of nodes using 
//' \code{\link{set_allele_freqs}}. Note that the initial allele distribution can also
//' be set in the next step, allowing this step to be skipped.
//' \item Simulate spread of genetic material through the network with 
//' \code{\link{popgen_dirichlet}}.
//' \item Draw samples from the simulated population using 
//' \code{\link{draw_isolates}}.}
//'
//' A worked example is available in the 'overview' vignette.
//'
//' @examples
//' # 1) fluid model
//' # this model is pretty resilient to numerical issues
//' inp <- c(0L, 0L, 1L, 2L, 3L, 1L, 5L)
//' outp <- c(1L, 2L, 3L, 3L, 4L, 4L, 2L)
//' rates <- c(1, 1.5, 0.5, 0.1, 1, 0.1, 0.5)
//' edgelist <- data.frame(inp, outp, rates)
//' ext <- data.frame(c(0L, 5L), c(0.5, 0.5))
//' net <- popsnetwork(edgelist, ext, 0.1)
//'
//' # 2) units model
//' # for this model absolute numbers are relevant
//' # note that whether we use integer or string node ids is completely arbitrary 
//' 
//' \dontrun{
//' # this will produce an error in popsnetwork since node C's output is larger than
//' # its input (300 > 250)
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(150, 100, 300))
//' ext <- data.frame(node=c("A", "B"), rate=c(300, 100), input=c(1000, 1000))
//' net <- popsnetwork(el, ext, spread_model="units")
//' }
//' 
//' # if we reduce C's output a bit it works fine
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(150, 100, 200))
//' ext <- data.frame(node=c("A", "B"), rate=c(300, 100), input=c(1000, 1000))
//' net <- popsnetwork(el, ext, spread_model="units")
//'
//' @param links A dataframe describing all edges in the graph as well as transfer rates
//' between them. The first two columns are read as inputs and outputs. If there are only 
//' two columns all rates are assumed to be 1.
//' @param external A dataframe describing external inputs into the network. The first column
//' is expected to contain node ids (as indices or factors), the second column specifies 
//' the (absolute) amount of infected material in the input. If there is a third column 
//' present it can
//' be used to set overall input rates on the respective nodes (this is relevant for the ibm).
//' If no rates are given an input of 1 is assumed for all input nodes. 
//' @param transmission Rate of infection within nodes (i.e. proportion of uninfected material
//' becoming infected).
//' @param decay The decay of material within nodes.
//' If this parameter has a value in [0, 1) transport rates for the entire network will 
//' be rescaled so that sum(output) == sum(input) * (1-decay) in all (non-leaf) nodes.
//' @param spread_model How to model spread of pathogens. With "fluid" the substrate carrying
//' the infection and the pathogen itself is essentially treated as a fluid and rates are 
//' calculated deterministically. With "units" infection as well as selection of infected vs.
//' uninfected material at outputs is modelled as a stochastic process on discrete units.
//' @param checks Perform some basic integrity checks on input data (currently looks for cycles
//' and disconnected sub-networks).
//' @return A popsnetwork object.
// [[Rcpp::export]]
XPtr<Net_t> popsnetwork(const DataFrame & links, const DataFrame & external, double transmission=0.0, double decay=-1.0, const string & spread_model = "fluid", bool checks=false);


// [[Rcpp::export(name=".printpopsnetwork")]]
void print_popsnetwork(const XPtr<Net_t> & p_net);


//' @title set_allele_freqs
//' 
//' @description Pre-set allele frequencies for some nodes.
//' 
//' @details Use this function to initialize allele frequencies for some nodes of the network.
//' 
//' @param p_net A popsnetwork object.
//' @param ini_dist Initial distribution of allele frequencies. ini_dist has to be 
//' a list
//' containing a vector of node IDs (see \code{\link{popsnetwork}}) and
//' a matrix of allele frequencies. Note that if non-root nodes are initialized in this way
//' they will be marked as blocked, i.e. they will not receive input in genetic material 
//' (but will produce output themselves).
//' @return A new popsnetwork object.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' set_allele_freqs(net, list(as.factor(c("A", "C")), freqs))
// [[Rcpp::export]]
XPtr<Net_t> set_allele_freqs(const XPtr<Net_t> & p_net, const List & ini_dist);


//' @title popgen_dirichlet
//' 
//' @description Simulate spread of pathogens on the network using a Dirichlet
//' distribution to approximate genetic drift.
//' 
//' @details This function simulates the change of gene frequencies in a population
//' of pathogens as they spread through the transport network starting at the 
//' external sources (see \code{\link{popsnetwork}}). At each node founder
//' effects are assumed to change composition of the population. This change is 
//' approximated by drawing a set of allele frequencies from a Dirichlet distribution.
//' 
//' @param p_net A popsnetwork object.
//' @param theta Scale parameter of the Dirichlet distribution. At each node the Dirichlet
//' distribution the new allele frequencies ar drawn from is parameterized by the old
//' frequencies multiplied by theta.
//' @param ini_dist Initial distribution of allele frequencies (optional). ini_dist has to be 
//' a list
//' containing a vector of node IDs (see \code{\link{popsnetwork}}) and
//' a matrix of allele frequencies. Note that *any* node pre-set in this
//' way will effectively be treated as a source and hide nodes that are further upstream (see
//' \code{\link{set_allele_freqs}}).
//' @return A new popsnetwork object with allele frequencies set for each node.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # run simulation on initialized net
//' ini_net <- set_allele_freqs(net, ini_freqs)
//' popgen_dirichlet(ini_net, 0.3)
//'
//' # or we can initialize and run in one call
//' popgen_dirichlet(net, 0.3, ini_freqs)
// [[Rcpp::export]]
XPtr<Net_t> popgen_dirichlet(const XPtr<Net_t> & p_net, double theta, Nullable<List> ini_dist = R_NilValue);


//' @title popgen_ibm_mixed
//' 
//' @description Simulate spread of pathogens on the network using a (very) simple individual-based
//' model.
//' 
//' @details This function simulates the change of gene frequencies in a population
//' of pathogens as they spread through the transport network starting at the 
//' external sources (see \code{\link{popsnetwork}}). At each node founder
//' effects are assumed to change composition of the population. This change is 
//' simulated by directly drawing from the distribution of genotypes and unfected units,
//' respectively. Note that for numerical reasons this model will not work for small
//' absolute numbers of units (input and transfer rates). It is also strongly recommended
//' to only run this model on networks that have been generated with the "units" method
//' (see \code{\link{popsnetwork}}).
//' 
//' @param p_net A popsnetwork object.
//' @param ini_dist Initial distribution of allele frequencies (optional). ini_dist has to be 
//' a list
//' containing a vector of node IDs (see \code{\link{popsnetwork}}) and
//' a matrix of allele frequencies. Note that *any* node pre-set in this
//' way will effectively be treated as a source and hide nodes that are further upstream (see
//' \code{\link{set_allele_freqs}}).
//' @return A new popsnetwork object with allele frequencies set for each node.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(150, 100, 200))
//' ext <- data.frame(node=c("A", "B"), rate=c(300, 100), input=c(1000, 1000))
//' net <- popsnetwork(el, ext, spread_model="units")
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # run simulation on initialized net
//' ini_net <- set_allele_freqs(net, ini_freqs)
//' popgen_ibm_mixed(ini_net)
//'
//' # or we can initialize and run in one call
//' popgen_ibm_mixed(net, ini_freqs)
// [[Rcpp::export]]
XPtr<Net_t> popgen_ibm_mixed(const XPtr<Net_t> & p_net, Nullable<List> ini_dist = R_NilValue);


//' @title draw_isolates
//'
//' @description Draw a set of isolates from the network.
//' 
//' @details Draw a random set of isolates from a number of nodes in the network. This will
//' *only* work if allele frequencies have been set or simulated.
//'
//' @param p_net a PopsNet object.
//' @param samples Number of samples to draw from each node. This has to be a dataframe
//' with node ids (see \code{\link{popsnetwork}}) in the first and number of isolates to 
//' draw in the second column.
//' @param aggregate Whether to return one line per sample taken or to sum up allele counts 
//' per node.
//' @return A dataframe containing 
//' \itemize{ 
//' \item if aggregate==TRUE: node id in \code{$node} and number of isolates with allele 
//' 	\code{x} in \code{$allele_x}.
//' \item if aggregate==FALSE: node id in \code{$node} and allele index in \code{$allele}}
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # simulate
//' res <- popgen_dirichlet(net, 0.3, ini_freqs)
//'
//' # get some data
//' draw_isolates(res, data.frame(c("C", "D"), c(10, 10)))
// [[Rcpp::export]]
DataFrame draw_isolates(const XPtr<Net_t> & p_net, const DataFrame & samples, 
	bool aggregate=true);


//' @title draw_alleles
//'
//' @description Draw a set of alleles from the network.
//' 
//' @details Draw a random set of alleles from a number of nodes in the network. This will
//' *only* work if allele frequencies have been set or simulated.
//'
//' @param p_net a PopsNet object.
//' @param nodes A vector of node ids (either integer or factor).
//' @param n How many alleles to draw per node.
//' @return A dataframe with one column per node containing a list of allele ids.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # simulate
//' res <- popgen_dirichlet(net, 0.3, ini_freqs)
//'
//' # get some data
//' draw_alleles(res, as.factor(c("C", "D")))
// [[Rcpp::export]]
DataFrame draw_alleles(const XPtr<Net_t> & p_net, const IntegerVector & nodes, int n=1);


//' @title egdeList
//'
//' @description Get a list of edges in a dataframe.
//'
//' @details This function returns a list of the edges in the network and the corresponding
//' transport rates.
//'
//' @param p_net A PopsNet object.
//' @param as_string Whether to return the list of nodes as string vector. If this is FALSE
//' (the default) a factor or plain integer vector (depending on how the net was constructed)
//' will be returned.
//' @return A dataframe with from, to, rates and rates_infected.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//' # get nodes as factor
//' edge_list(net)
//' # get nodes as string
//' edge_list(net, TRUE)
// [[Rcpp::export]]
DataFrame edge_list(const XPtr<Net_t> & p_net, bool as_string=false);

//' @title node_list
//'
//' @description Get a list of nodes in a dataframe.
//'
//' @details This function returns a list of the nodes in the network as well as the amount of
//' infected material they contain.
//'
//' @param p_net A PopsNet object.
//' @param as_string Whether to return the list of nodes as string vector. If this is FALSE
//' (the default) a factor or plain integer vector (depending on how the net was constructed)
//' will be returned.
//' @return A dataframe with id and rate_infected.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//' node_list(net)
// [[Rcpp::export]]
DataFrame node_list(const XPtr<Net_t> & p_net, bool as_string=false);

//' @title distances_topology
//'
//' @description Calculate topological distances between nodes in a network.
//' 
//' @details This function calculates the topological distance (number of edges in
//' the shortest path) between all pairs of nodes in a network.
//' 
//' @param p_net A popsnetwork object.
//' @param leaves_only Whether to save time by generating only distances between leave nodes.
//' The rest of the distance matrix will be filled with -1 or 0 (diagonal) in this case.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # get distances
//' distances_topology(net)
// [[Rcpp::export]]
NumericMatrix distances_topology(const XPtr<Net_t> & p_net, bool leaves_only = true);

//' @title distances_freqdist
//'
//' @description Calculate genetic dissimilarities within a network.
//' 
//' @details This function calculates the dissimilarity (mean square distance in
//' allele frequencies) of all pairs of nodes in a network.
//' 
//' @param p_net A popsnetwork object.
//' @param skip_empty Whether to return NA for empty nodes.
//' @return A matrix of all distances.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # simulate
//' res <- popgen_dirichlet(net, 0.3, ini_freqs)
//'
//' # get distances
//' distances_freqdist(res)
// [[Rcpp::export]]
NumericMatrix distances_freqdist(const XPtr<Net_t> & p_net, bool skip_empty=true);


//' @title distances_sample
//'
//' @description Calculate genetic distances within a network.
//' 
//' @details This function calculates the genetic distance of all pairs of nodes in a 
//' network by comparing a number of random samples from each node (using Hamming distance).
//' 
//' @param p_net A popsnetwork object.
//' @param n How many samples per node to use for comparison.
//' @param skip_empty Whether to return NA for empty nodes.
//' @return A matrix of all distances.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # simulate
//' res <- popgen_dirichlet(net, 0.3, ini_freqs)
//'
//' # get distances
//' distances_sample(res)
// [[Rcpp::export]]
NumericMatrix distances_sample(const XPtr<Net_t> & p_net, int n=1, bool skip_empty=true);


//' @title distances_EHamming
//'
//' @description Calculate genetic distances within a network as expected values of Hamming
//' distances.
//' 
//' @details This function calculates the genetic distance of all pairs of nodes in a 
//' network by calculating per pair of nodes the average Hamming distance between them 
//' (more precisely the expected value of the Hamming distance between two individuals 
//' randomly selected from each of the nodes).
//' 
//' @param p_net A popsnetwork object.
//' @param skip_empty Whether to return NA for empty nodes.
//' @return A matrix of all distances.
//'
//' @examples
//' # create network
//' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
//' ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))
//' net <- popsnetwork(el, ext)
//'
//' # set allele frequencies (2 nodes, 3 alleles)
//' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
//' ini_freqs <- list(as.factor(c("A", "C")), freqs)
//'
//' # simulate
//' res <- popgen_dirichlet(net, 0.3, ini_freqs)
//'
//' # get distances
//' distances_EHamming(res)
// [[Rcpp::export]]
NumericMatrix distances_EHamming(const XPtr<Net_t> & p_net, bool skip_empty=true);


//' @title generate_PA
//'
//' @description Generate a random transport network using preferential attachment.
//'
//' @details This function generates a random scale-free network. It uses the Barabasi-Albert 
//' preferential attachment algorithm, slightly modified to allow for directedness and 
//' isolated initial nodes.
//'
//' @param n_nodes Number of (non-source) nodes to generate.
//' @param n_sources Number of source nodes to initialize the network with (has to be at least
//' 1). Note that there is no guarantee all source nodes will become part of the network.
//' @param m_dist The probability distribution to draw the number of inputs for new nodes
//' from. m_dist will be normalized, therefore it does not have to sum up to 1.
//' @param zero_appeal Constant to be added to the nodes' attractiveness.
//' @param compact Whether to remove isolated source nodes.
//' @return An edgelist as a dataframe. 
//'
//' @examples
//' generate_PA(20, 5, c(3, 1))
// [[Rcpp::export]]
DataFrame generate_PA(int n_nodes, int n_sources, NumericVector m_dist, float zero_appeal=1,
	bool compact=true);


#endif	// DIR_NETWORK_H
