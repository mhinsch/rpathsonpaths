#include <vector>
#include <iostream>

#include "genericgraph.h"
#include "transportgraph.h"
#include "driftapprox.h"
#include "network.h"
#include "network_io.h"

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
	vector<double> scaled;
	gsl_rng * rng;
	double theta;

	Drift(size_t size, double t)
		: theta(t)
		{
		rng = gsl_rng_alloc(gsl_rng_default);
		}

	void operator()(const N_t::freq_t & freqs, N_t::freq_t & result)
		{
		assert(result.size() == freqs.size());
		scaled.resize(freqs.size());

		size_t i=0;
		for (auto f : freqs)
			scaled[i++] = f * theta;

		gsl_ran_dirichlet(rng, scaled.size(), scaled.data(), result.data());
		}
	};


int main()
	{
	Net_t net;
	read_network(cin, net);
	
	gsl_rng_env_setup();

	Drift drift(4, 100000);

	size_t i = 0;
	for (auto & n : net.nodes)
		{
		assert(n);
		assert(n->consistent());
		if (n->is_root())
			{
			n->frequencies.resize(4, 0);
			n->frequencies[i++] = 1.0;
			}
		}

	annotate_rates(net.nodes.begin(), net.nodes.end(), 0.01);
	annotate_frequencies(net.nodes.begin(), net.nodes.end(), drift);

	i = 0;
	for (auto n : net.nodes)
		{
		cout << i++ << ":\t" << n->rate_in << "\t" << n->rate_in_infd << "\n";
		
		for (auto f : n->frequencies)
			cout << "\t" << f;
		cout << "\n";
		//assert(n->valid(0.0001));
		}

	cout << "cloning...\n";

	i = 0;

	Net_t * second = net.clone();
	for (auto n : second->nodes)
		{
		cout << i++ << ":\t" << n->rate_in << "\t" << n->rate_in_infd << "\n";
		
		for (auto f : n->frequencies)
			cout << "\t" << f;
		cout << "\n";
		//assert(n->valid(0.0001));
		}

	delete second;
	}
