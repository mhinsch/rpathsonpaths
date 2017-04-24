#ifndef DIR_NETWORK_H
#define DIR_NETWORK_H

#include "libpathsonpaths/genericgraph.h"
#include "libpathsonpaths/transportgraph.h"
#include "libpathsonpaths/driftapprox.h"
#include "libpathsonpaths/network.h"

#include <vector>

#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

template<class T>
using StdVector = vector<T>;

template<class GRAPH>
using MyVecDriftNode = DriftNode<vector<double>, TranspNode<Node<GRAPH, StdVector> > >;

template<class GRAPH>
using MyTranspLink = TranspLink<Link<GRAPH> >;

typedef Graph<MyVecDriftNode, MyTranspLink> G_t;
typedef G_t::node_t N_t;
typedef G_t::link_t L_t;

typedef Network<N_t, L_t> Net_t;

struct Drift
	{
	typedef typename N_t::freq_t::data_type num_t;
	vector<num_t> result;
	num_t theta;

	Drift(size_t size, double t)
		: result(size), theta(t)
		{ }

	void operator()(N_t::freq_t & freqs)
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

// [[Rcpp::export]]
XPtr<Net_t> PopsNetwork(List links, List external, double transmission);

// [[Rcpp::export,name=(".printPopsNetwork"]]
void print_PopsNetwork(const XPtr<Net_t> & pNet);

// [[Rcpp::export,name=(".printPopsNode"]]
void print_PopsNode(const XPtr<N_t> & pNode);

// [[Rcpp::export]]
XPtr<Net_t> spreadDirichlet(const XPtr<Net_t> & pNet, List iniDist, double theta);

// [[Rcpp::export(name="getPopsNode")]]
XPtr<N_t> getPopsNode(const XPtr<Net_t> & pNet, int id);

// [[Rcpp::export(name="drawIsolates.PopsNode")]]
List drawIsolates_PopsNode(const XPtr<N_t> & pNode, int n);

// [[Rcpp::export(name="drawIsolates.PopsNetwork")]]
List drawIsolates_PopsNetwork(const XPtr<N_t> & pNet, List samples);

#endif	// DIR_NETWORK_H
