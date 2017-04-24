#include "libpathsonpaths/genericgraph.h"
#include "libpathsonpaths/transportgraph.h"
#include "libpathsonpaths/driftapprox.h"
#include "libpathsonpaths/network.h"
#include "libpathsonpaths/network_io.h"

#include <vector>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

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

// TODO use R rng
struct Drift
	{
	vector<double> result;
	vector<double> scaled;
	gsl_rng * rng;
	double theta;

	Drift(size_t size, double t)
		: result(size), theta(t)
		{
		rng = gsl_rng_alloc(gsl_rng_default);
		}

	void operator()(const N_t::freq_t & freqs)
		{
		assert(result.size() == freqs.size());
		scaled.resize(freqs.size());

		size_t i=0;
		for (auto f : freqs)
			scaled[i++] = f * theta;

		gsl_ran_dirichlet(rng, scaled.size(), scaled.data(), result.data());
		}

	vector<double>::const_iterator begin() const
		{
		return result.begin();
		}
	
	vector<double>::const_iterator end() const
		{
		return result.end();
		}
	};

// [[Rcpp::export]]
XPtr<Net_t> PopsNetwork(List links, List external, double transmission)
	{
	const IntegerVector inputs = links["inputs"];
	const NumericVector rates = links["rates"];
	const IntegerVector outputs = links["outputs"];

	// TODO check sizes of arrays

	IntegerVector ext = external["nodes"];
	NumericVector ext_rates = external["rates"];

	// TODO check sizes of arrays
	

	Net_t * net = new Net_t;

	for (size_t i=0; i<inputs.size(); i++)
		net->add_link(inputs[i], rates[i], outputs[i]);

	for (size_t i=0; i<external.size(); i++)
		net->set_source(external[i], ext_rates[i]);

	annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);

	return XPtr<Net_t>(net, true);
	}

// [[Rcpp::export,name=("format.PopsNetwork"]]
List format_PopsNetwork(const XPtr<Net_t> & pNet)
	{
	Net_t * net = pNet.get();
	}

// [[Rcpp::export,name=("format.PopsNode"]]
List format_PopsNode(const XPtr<N_t> & pNode)
	{
	N_t * node = pNode.get();
	}

// [[Rcpp::export]]
XPtr<Net_t> spreadDirichlet(const XPtr<Net_t> & pNet, List iniDist, double theta)
	{
	Net_t * net = p_net.get().clone();
	const IntegerVector nodes = iniDist["nodes"];
	const NumericMatrix freqs = IniDist["frequencies"];

	const size_t n_all = freqs.ncol();

	// TODO check sizes of arrays
	
	for (size_t i=0; i<nodes.size(); i++)
		{
		net->nodes[nodes[i]].frequencies.resize(n_all);
		for (size_t j=0; j<n_all; j++)
			net->nodes[nodes[i]].frequencies[j] = freqs(i, j);
		}

	// TODO rng

	Drift drift(n_all, theta);
	annotate_frequencies(net.nodes.begin(), net.nodes.end(), drift);
	
	return XPtr<Net_t>(net, true);
	}

// [[Rcpp::export(name="getNode.PopsNetwork")]]
XPtr<N_t> getNode_PopsNetwork(const XPtr<Net_t> & pNet, int id)
	{
	Net_t * net = pNet.get();

	// TODO check validity of id
	
	N_t * node = net->nodes[i];

	// don't GC, since net owns the memory
	return XPtr<N_t>(node, false);
	}

// [[Rcpp::export(name="drawIsolates.PopsNode")]]
List drawIsolates_PopsNode(const XPtr<N_t> & pNode, int n)
	{
	N_t * node = pNode.get();
	}

// [[Rcpp::export(name="drawIsolates.PopsNetwork")]]
List drawIsolates_PopsNetwork(const XPtr<N_t> & pNet, List samples)
	{
	Network * net = pNet.get();

	IntegerVector nodes = samples["nodes"];
	IntegerVector num = samples["N"];

	// TODO check array sizes	
	
	}

