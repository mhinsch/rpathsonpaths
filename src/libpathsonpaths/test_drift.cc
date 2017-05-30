#include <vector>
#include <iostream>

#include "genericgraph.h"
#include "transportgraph.h"
#include "driftapprox.h"
#include "transportnetwork.h"
#include "network_io.h"
#include "ibmmixed.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



using namespace std;

template<class T>
using StdVector = vector<T>;


template<class GRAPH>
struct MyDriftNode : 
	public FreqNode<vector<double>>, 
	public TranspNode,
	public Node<GRAPH, StdVector>
	{};

template<class GRAPH>
struct MyTranspLink : public TranspLink, public Link<GRAPH>
	{
	MyTranspLink(typename GRAPH::node_t * f, typename GRAPH::node_t * t,
		double a_rate = 0.0, double a_rate_infd = -1)
		: TranspLink(a_rate, a_rate_infd), Link<GRAPH>(f, t)
		{}
	};


typedef Graph<MyDriftNode, MyTranspLink> G_t;
typedef G_t::node_t N_t;
typedef G_t::link_t L_t;
typedef TransportNetwork<N_t, L_t> Net_t;


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


struct Rng
	{
	gsl_rng * rng;

	Rng()
		{
		rng = gsl_rng_alloc(gsl_rng_default);
		}

	int binom(double p, int t)
		{
		return gsl_ran_binomial(rng, p, t);
		}

	int hypergeom(double n1, double n2, double k)
		{
		return gsl_ran_hypergeometric(rng, n1, n2, k);
		}
	};


int main()
	{
	Net_t net;
	read_network(cin, net);
	
	gsl_rng_env_setup();

	Drift drift(4, 10);

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

	Net_t net2 = net;

	annotate_frequencies(net.nodes.begin(), net.nodes.end(), drift);

	i = 0;
	for (auto n : net.nodes)
		{
		cout << i++ << ": " << n->rate_in << ", " << n->rate_in_infd << ", "
			<< n->d_rate_in_infd << ", "
			<< n->rate_out_infd << "\n";
		
		for (auto f : n->frequencies)
			cout << "\t" << f;
		cout << "\n\n";
		//assert(n->valid(0.0001));
		}

	cout << "ibm...\n";

	Rng rng;

	for (auto & n : net2.nodes)
		{
		if (n->is_root())
			{
			n->rate_in = 1000;
			for (auto & f : n->frequencies)
				f *= 100;
			}
		}
	
	annotate_frequencies_ibmm(net2.nodes.begin(), net2.nodes.end(), rng);

	i = 0;
	for (auto n : net2.nodes)
		{
		cout << i++ << ": " << n->rate_in << ", " << n->rate_in_infd << ", "
			<< n->d_rate_in_infd << ", "
			<< n->rate_out_infd << "\n";
		
		for (auto f : n->frequencies)
			cout << "\t" << f;
		cout << "\n\n";
		//assert(n->valid(0.0001));
		}
	}
