#ifndef DIR_NETWORK_H
#define DIR_NETWORK_H

#include "rpathsonpaths_types.h"

struct Drift
	{
	typedef typename Node_t::freq_t::value_type num_t;
	vector<num_t> result;
	num_t theta;

	Drift(size_t size, double t)
		: result(size), theta(t)
		{ }

	void operator()(Node_t::freq_t & freqs)
		{
		assert(result.size() == freqs.size());
		num_t norm = 0.0;		

		for (num_t & f : freqs)
			norm += (f = R::rgamma(f * theta, 1.0));

		for (num_t & f : freqs)
			f /= norm;
		}

	vector<num_t>::const_iterator begin() const
		{
		return result.begin();
		}
	
	vector<num_t>::const_iterator end() const
		{
		return result.end();
		}
	};

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
XPtr<Net_t> PopsNetwork(DataFrame links, DataFrame external, double transmission);

// [[Rcpp::export,name=(".printPopsNetwork"]]
void print_PopsNetwork(const XPtr<Net_t> & pNet);

// [[Rcpp::export,name=(".printPopsNode"]]
void print_PopsNode(const XPtr<Node_t> & pNode);

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
//' @param iniDist Initial distribution of allele frequencies. iniDist has to be a list
//' containing a vector of node IDs (@seealso \code{\link{PopsNetwork}}) as $nodes and
//' a matrix of allele frequencies as $frequencies. Note that *any* node pre-set in this
//' way will be treated as a source only.
//' @param theta Scale parameter of the Dirichlet distribution. At each node the Dirichlet
//' distribution the new allele frequencies ar drawn from is parameterized by the old
//' frequencies multiplied by theta.
//' @return A new PopsNetwork object with allele frequencies set for each node.
// [[Rcpp::export]]
XPtr<Net_t> spreadDirichlet(const XPtr<Net_t> & pNet, const List iniDist, double theta);

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
DataFrame drawIsolates_PopsNetwork(const XPtr<Net_t> & pNet, DataFrame samples);

#endif	// DIR_NETWORK_H
