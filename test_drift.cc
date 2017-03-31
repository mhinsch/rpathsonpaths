#include "genericgraph.h"
#include "transportgraph.h"
#include "driftapprox.h"
#include "network.h"
#include "network_io.h"

#include <vector>
#include <iostream>

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

int main()
	{
	Net_t net;
	read_network(cin, net);
	
	annotate_rates(net.nodes.begin(), net.nodes.end(), 0.01);
	}
