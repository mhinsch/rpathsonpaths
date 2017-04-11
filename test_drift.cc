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


int main()
	{
	Net_t net;
	read_network(cin, net);
	
	gsl_rng_env_setup();

	Drift drift(4, 10);

	for (auto & n : net.nodes)
		{
		assert(n);
		assert(n->consistent());
		if (n->is_root())
			n->frequencies.resize(4, 0.25);
		}

	annotate_rates(net.nodes.begin(), net.nodes.end(), 0.01);
	annotate_frequencies(net.nodes.begin(), net.nodes.end(), drift);

	size_t i = 0;
	for (auto n : net.nodes)
		{
		assert(n->valid(0.0001));
		cout << i << "\t" << n->rate_in << "\t" << n->rate_in_infd << "\n";
		for (auto f : n->frequencies)
			cout << "\t" << f;
		cout << "\n";
		}
	}
