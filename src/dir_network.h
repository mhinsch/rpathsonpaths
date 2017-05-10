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
//' \code{\link{spread_dirichlet}}.
//' \item Draw samples from the simulated population using \code{\link{draw_isolates}}.}
//'
//' A worked example is available in the 'overview' vignette.
//'
//' @examples
//'
//' inp <- c(0L, 0L, 1L, 2L, 3L, 1L, 5L)
//' outp <- c(1L, 2L, 3L, 3L, 4L, 4L, 2L)
//' rates <- c(1, 1.5, 0.5, 0.1, 1, 0.1, 0.5)
//' edgelist <- data.frame(inp, outp, rates)
//' ext <- data.frame(c(0L, 5L), c(0.5, 0.5))
//' net <- popsnetwork(edgelist, ext, 0.1)
//' 
//'
//' @param links A dataframe describing all edges in the graph as well as transfer rates
//' between them. The first two columns are read as inputs and outputs. If there are only 
//' two columns all rates are assumed to be 1.
//' @param external A dataframe describing external inputs into the network. The first column
//' is expected to contain node ids (as indices or factors), the second column specifies 
//' the proportion of infected material in the input.
//' @param transmission Rate of infection within nodes (i.e. proportion of uninfected material
//' becoming infected).
//' @param checks Perform some basic integrity checks on input data (currently looks for cycles
//' and disconnected sub-networks).
//' @return A popsnetwork object.
// [[Rcpp::export]]
XPtr<Net_t> popsnetwork(const DataFrame & links, const DataFrame & external, double transmission=0.0, bool checks=false);

// [[Rcpp::export(name=".printpopsnetwork")]]
void print_popsnetwork(const XPtr<Net_t> & p_net);

// [[Rcpp::export(name=".printpopsnode")]]
void print_popsnode(const XPtr<Node_t> & p_node);

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
//' a matrix of allele frequencies. Note that *any* node pre-set in this
//' way will be treated as a source by \code{spread_dirichlet} (this hiding nodes further
//' upstream).
//' @return A new popsnetwork object.
// [[Rcpp::export]]
XPtr<Net_t> set_allele_freqs(const XPtr<Net_t> & p_net, const List & ini_dist);


//' @title spread_dirichlet
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
//' way will effectively be treated as a source and hide nodes that are further upstream.
//' @return A new popsnetwork object with allele frequencies set for each node.
// [[Rcpp::export]]
XPtr<Net_t> spread_dirichlet(const XPtr<Net_t> & p_net, double theta, Nullable<List> ini_dist = R_NilValue);

//' @title get_popsnode
//'
//' @description Pick a single node from the network.
//' 
//' @param p_net A PopsNet object.
//' @param id The id of the node to return (either a string or an integer, see 
//' \code{\link{popsnetwork}}).
//' @return A popsnode object.
// [[Rcpp::export(name="get_popsnode")]]
XPtr<Node_t> get_popsnode(const XPtr<Net_t> & p_net, SEXP id);

//' @title draw_isolates.popsnode
//'
//' @description Draw a set of isolates from a single node.
//'
//' @details Draw a random set of isolates from a given node. This will only work if
//' allele frequencies have been set (manually or by simulation).
//'
//' @param p_node A popsnode object.
//' @param n Number of isolates to draw.
//' @return A vector of allele counts.
// [[Rcpp::export(name="draw_isolates.popsnode")]]
IntegerVector draw_isolates_popsnode(const XPtr<Node_t> & p_node, int n);

//' @title draw_isolates.popsnetwork
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
//' @return A dataframe containing node id in $node and number of isolates with allele \code{x}
//' in $\code{allele_x}.
// [[Rcpp::export(name="draw_isolates.popsnetwork")]]
DataFrame draw_isolates_popsnetwork(const XPtr<Net_t> & p_net, const DataFrame & samples);

//' @title egdeList
//'
//' @description Get a list of edges in a dataframe.
//'
//' @description This function returns a list of the edges in the network in a format
//' that is suitable for plotting with e.g. iGraph.
//'
//' @param p_net A PopsNet object.
//' @return A dataframe with from, to, rates and rates_infected.
// [[Rcpp::export]]
DataFrame edge_list(const XPtr<Net_t> & p_net);

//' @title node_list
//'
//' @description Get a list of nodes in a dataframe.
//'
//' @description This function returns a list of the nodes in the network in a format
//' that is suitable for plotting with e.g. iGraph.
//'
//' @param p_net A PopsNet object.
//' @return A dataframe with id and rate_infected.
// [[Rcpp::export]]
DataFrame node_list(const XPtr<Net_t> & p_net);

#endif	// DIR_NETWORK_H
