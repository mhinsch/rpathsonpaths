#ifndef DIR_NETWORK_H
#define DIR_NETWORK_H

#include "rpathsonpaths_types.h"

struct Drift
	{
	typedef typename Node_t::freq_t::value_type num_t;
	num_t theta;

	Drift(double t)
		: theta(t)
		{ }

	void operator()(const Node_t::freq_t & freqs, Node_t::freq_t & res)
		{
		assert(res.size() == freqs.size());
		num_t norm = 0.0;		

		for (size_t i=0; i<freqs.size(); i++)
			norm += (res[i] = R::rgamma(freqs[i] * theta, 1.0));

		for (num_t & f : res)
			f /= norm;
		}
	};



//' @title sources
//'
//' @description Get a list of source nodes from a network.
//'
//' @details Extract a list of source nodes, i.e. nodes that do not have an 
//' input, from a network.
//'
//' @param edgeList A dataframe containing a list of edges.
//'
//' @return An integer vector containg the ids of all source nodes in the network.
// [[Rcpp::export]]
IntegerVector sources(const DataFrame & edgeList);

//' @title colourNetwork
//'
//' @description Identify separate sub-networks.
//'
//' @details This function identifies completely separate sub-networks in a network
//' described as an edge list.
//'
//' @param edgeList A dataframe containing a list of edges.
//'
//' @return An integer vector with the sub-network id of each edge. Note that id's start at
//' 1 and are not guaranteed to be contiguous.
// [[Rcpp::export]]
IntegerVector colourNetwork(const DataFrame & edgeList);
	
//' @title PopsNetwork 
//'
//' @description Create a PopsNetwork object.
//'
//' @details A PopsNetwork object stores the nodes and edges making up a food transport
//' network and associated data describing spread of infected diseases on the
//' network. A network is created from a tabular description of its edges.
//'
//' @param links A dataframe describing all edges in the graph in the columns 
//' $inputs, $outputs and $rates. Note that inputs and outputs refer to node ids.
//' @param external A dataframe describing external inputs into the network as $nodes
//' and $rates.
//' @param transmission Rate of infection within nodes.
//' @return A PopsNetwork object.
// [[Rcpp::export]]
XPtr<Net_t> PopsNetwork(const DataFrame & links, const DataFrame & external, double transmission);

// [[Rcpp::export(name=".printPopsNetwork")]]
void print_PopsNetwork(const XPtr<Net_t> & pNet);

// [[Rcpp::export(name=".printPopsNode")]]
void print_PopsNode(const XPtr<Node_t> & pNode);

//' @title setAlleleFreqs
//' 
//' @description Pre-set allele frequencies for some nodes.
//' 
//' @details Use this function to initialize allele frequencies for some nodes of the network.
//' 
//' @param pNet A PopsNetwork object.
//' @param iniDist Initial distribution of allele frequencies. iniDist has to be 
//' a list
//' containing a vector of node IDs (@seealso \code{\link{PopsNetwork}}) as $nodes and
//' a matrix of allele frequencies as $frequencies. Note that *any* node pre-set in this
//' way will be treated as a source only by spreadDirichlet.
//' @return A new PopsNetwork object.
// [[Rcpp::export]]
XPtr<Net_t> setAlleleFreqs(const XPtr<Net_t> & pNet, const List & iniDist);


//' @title spreadDirichlet
//' 
//' @description Simulate spread of pathogens on the network using a Dirichlet
//' distribution to approximate genetic drift.
//' 
//' @details This function simulates the change of gene frequencies in a population
//' of pathogens as they spread through the transport network starting at the 
//' external sources (@seealso \code{\link{PopsNetwork}}). At each node founder
//' effects are assumed to change composition of the population. This change is 
//' approximated by drawing a set of allele frequencies from a Dirichlet distribution.
//' 
//' @param pNet A PopsNetwork object.
//' @param theta Scale parameter of the Dirichlet distribution. At each node the Dirichlet
//' distribution the new allele frequencies ar drawn from is parameterized by the old
//' frequencies multiplied by theta.
//' @param iniDist Initial distribution of allele frequencies (optional). iniDist has to be 
//' a list
//' containing a vector of node IDs (@seealso \code{\link{PopsNetwork}}) as $nodes and
//' a matrix of allele frequencies as $frequencies. Note that *any* node pre-set in this
//' way will be treated as a source only.
//' @return A new PopsNetwork object with allele frequencies set for each node.
// [[Rcpp::export]]
XPtr<Net_t> spreadDirichlet(const XPtr<Net_t> & pNet, double theta, Nullable<List> iniDist = R_NilValue);

//' @title getPopsNode
//'
//' @description Pick a single node from the network.
//' 
//' @param pNet A PopsNet object.
//' @param id The id of the node to return.
//' @return A PopsNode object.
// [[Rcpp::export(name="getPopsNode")]]
XPtr<Node_t> getPopsNode(const XPtr<Net_t> & pNet, int id);

//' @title drawIsolates.PopsNode
//'
//' @description Draw a set of isolates from a single node.
//'
//' @details Draw a random set of isolates from a given node. This will only work if
//' allele frequencies have been set (manually or by simulation).
//'
//' @param pNode A PopsNode object.
//' @param n Number of isolates to draw.
//' @return A vector of allele counts.
// [[Rcpp::export(name="drawIsolates.PopsNode")]]
IntegerVector drawIsolates_PopsNode(const XPtr<Node_t> & pNode, int n);

//' @title drawIsolates.PopsNetwork
//'
//' @description Draw a set of isolates from the network.
//' 
//' @details Draw a random set of isolates from a number of nodes in the network. This will
//' *only* work if allele frequencies have been set or simulated.
//'
//' @param pNet a PopsNet object.
//' @param samples Number of samples to draw from each node. This has to be a dataframe
//' with node ids in $nodes and number of isolates to draw in $N.
//' @return A dataframe with node id in $node and number of isolates with allele \code{x}
//' in $\code{allele_x}.
// [[Rcpp::export(name="drawIsolates.PopsNetwork")]]
DataFrame drawIsolates_PopsNetwork(const XPtr<Net_t> & pNet, const DataFrame & samples);

//' @title egdeList
//'
//' @description Get a list of edges in a dataframe.
//'
//' @description This function returns a list of the edges in the network in a format
//' that is suitable for plotting with e.g. iGraph.
//'
//' @param pNet A PopsNet object.
//' @return A dataframe with from, to, rates and rates_infected.
// [[Rcpp::export]]
DataFrame edgeList(const XPtr<Net_t> & pNet);

//' @title nodeList
//'
//' @description Get a list of nodes in a dataframe.
//'
//' @description This function returns a list of the nodes in the network in a format
//' that is suitable for plotting with e.g. iGraph.
//'
//' @param pNet A PopsNet object.
//' @return A dataframe with id and rate_infected.
// [[Rcpp::export]]
DataFrame nodeList(const XPtr<Net_t> & pNet);

#endif	// DIR_NETWORK_H
