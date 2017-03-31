#include "genericgraph.h"
#include "transportgraph.h"
#include "driftapprox.h"
#include "network.h"
#include "network_io.h"

#include <vector>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



using namespace std;

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
	vector<double> result;
	gsl_rng * rng;

	Drift(size_t size)
		: result(size)
		{
		rng = gsl_rng_alloc(gsl_rng_default);
		}

	void operator()(const N_t::freq_t & freqs)
		{
		assert(result.size() == freqs.size());
		gsl_ran_dirichlet(rng, freqs.size(), freqs.data(), result.data());
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


int main()
	{
	Net_t net;
	read_network(cin, net);
	
	gsl_rng_env_setup();

	Drift drift(10);

	for (auto & n : net.nodes)
		{
		assert(n);
		assert(n->consistent());
		if (n->is_root())
			n->frequencies.resize(10, 0.1);
		}

	annotate_rates(net.nodes.begin(), net.nodes.end(), 0.01);
	annotate_frequencies(net.nodes.begin(), net.nodes.end(), drift);
	}
