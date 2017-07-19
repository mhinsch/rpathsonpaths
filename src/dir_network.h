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
//' \code{\link{popgen_dirichlet}}.
//' \item Draw samples from the simulated population using 
//' \code{\link{draw_isolates.popsnetwork}}.}
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
//' the amount of infected material in the input. If there is a third column present it
//' be used to set overall input rates on the respective nodes (this is relevant for the ibm).
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
//' a matrix of allele frequencies. Note that if non-root nodes are initialized in this way
//' they will be marked as blocked, i.e. they will not receive input in genetic material 
//' (but will produce output themselves).
//' @return A new popsnetwork object.
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
//' way will effectively be treated as a source and hide nodes that are further upstream.
//' @return A new popsnetwork object with allele frequencies set for each node.
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
//' respectively.
//' 
//' @param p_net A popsnetwork object.
//' @param ini_dist Initial distribution of allele frequencies (optional). ini_dist has to be 
//' a list
//' containing a vector of node IDs (see \code{\link{popsnetwork}}) and
//' a matrix of allele frequencies. Note that *any* node pre-set in this
//' way will effectively be treated as a source and hide nodes that are further upstream.
//' @return A new popsnetwork object with allele frequencies set for each node.
// [[Rcpp::export]]
XPtr<Net_t> popgen_ibm_mixed(const XPtr<Net_t> & p_net, Nullable<List> ini_dist = R_NilValue);


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
//' @details This function returns a list of the edges in the network in a format
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
//' @details This function returns a list of the nodes in the network in a format
//' that is suitable for plotting with e.g. iGraph.
//'
//' @param p_net A PopsNet object.
//' @return A dataframe with id and rate_infected.
// [[Rcpp::export]]
DataFrame node_list(const XPtr<Net_t> & p_net);

//' @title SNP_distance
//'
//' @description Calculate distance (in number of different SNPs) between two genotypes 
//' encoded as integers. 
//' 
//' @details This function calculates the difference between two genetic sequences. A sequence
//' is represented as an integer number where the binary of that number corresponds
//' to a list of alleles. The distance between to sequences is the number of positions at
//' which they have different alleles.
//' 
//' @param g1 Sequence 1
//' @param g2 Sequence 2
//' @return The distance between g1 and g2.
// [[Rcpp::export]]
int SNP_distance(int g1, int g2);

//' @title SNP_distance_pop
//'
//' @description Calculate genetic distance between two populations.
//' 
//' @details This function calculates the genetic distance between two populations as the
//' average distance between all pairs of individuals of the two populations where genotype
//' is encoded as an integer.
//' 
//' @param p1 Population 1.
//' @param p2 Population 2.
//' @return The distance between p1 and p2.
// [[Rcpp::export]]
double SNP_distance_pop(const IntegerVector & p1, const IntegerVector & p2);

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
// [[Rcpp::export]]
NumericMatrix distances_EHamming(const XPtr<Net_t> & p_net, bool skip_empty=true);

#endif	// DIR_NETWORK_H
