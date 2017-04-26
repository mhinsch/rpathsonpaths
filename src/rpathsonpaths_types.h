#ifndef RPATHSONPATHS_TYPES_H
#define RPATHSONPATHS_TYPES_H

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
typedef G_t::node_t Node_t;
typedef G_t::link_t Link_t;

typedef Network<Node_t, Link_t> Net_t;


#endif	// RPATHSONPATHS_TYPES_H
