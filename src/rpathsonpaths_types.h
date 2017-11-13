#ifndef RPATHSONPATHS_TYPES_H
#define RPATHSONPATHS_TYPES_H


/** @file
 * Here we plug all the various template types from pathonpaths together to
 * create the types we need for the package. 
 */

#include "libpathsonpaths/genericgraph.h"
#include "libpathsonpaths/transportgraph.h"
#include "libpathsonpaths/driftapprox.h"
#include "libpathsonpaths/genefreqgraph.h"

#include "rnetwork.h"

#include <vector>


using namespace std;


// for some reason we need this
template<class T>
using StdVector = vector<T>;


// assemble all required node components
template<class GRAPH>
struct MyDriftNode : 
	public FreqNode<vector<double>>, 
	public TranspNode,
	public Node<GRAPH, StdVector>
	{};


// we need a new constructor, so we actually have to implement a new class here
template<class GRAPH>
struct MyTranspLink : public TranspLink, public Link<GRAPH>
	{
	MyTranspLink(typename GRAPH::node_t * f, typename GRAPH::node_t * t,
		double a_rate = 0.0, double a_rate_infd = 0)
		: TranspLink(a_rate, a_rate_infd), Link<GRAPH>(f, t)
		{}
	};


// tie everything together
typedef Graph<MyDriftNode, MyTranspLink> G_t;
typedef G_t::node_t Node_t;
typedef G_t::link_t Link_t;


// our network type
typedef RNetwork<Node_t, Link_t> Net_t;


#endif	// RPATHSONPATHS_TYPES_H
