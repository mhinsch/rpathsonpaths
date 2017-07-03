#ifndef RPATHSONPATHS_TYPES_H
#define RPATHSONPATHS_TYPES_H

#include "libpathsonpaths/genericgraph.h"
#include "libpathsonpaths/transportgraph.h"
#include "libpathsonpaths/driftapprox.h"
#include "libpathsonpaths/genefreqgraph.h"

#include "rnetwork.h"

#include <vector>


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
		double a_rate = 0.0, double a_rate_infd = 0)
		: TranspLink(a_rate, a_rate_infd), Link<GRAPH>(f, t)
		{}
	};


typedef Graph<MyDriftNode, MyTranspLink> G_t;
typedef G_t::node_t Node_t;
typedef G_t::link_t Link_t;

typedef RNetwork<Node_t, Link_t> Net_t;


#endif	// RPATHSONPATHS_TYPES_H
