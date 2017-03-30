#ifndef GENERICGRAPH_H
#define GENERICGRAPH_H

template<class GRAPH>
struct Link
	{
	typedef GRAPH graph_t;
	typedef graph_t::node_t node_t;

	node_t * from, * to;

	Link(node_t * f = 0, node_t * t = 0)
		: from(f), to(t)
		{}
	};

template<class GRAPH, template<class> class CONT>
struct Node
	{
	typedef GRAPH graph_t;
	typedef graph_t::link_t link_t;

	typedef CONT cont_t;

	cont_t<link_t *> inputs;
	cont_t<link_t *> outputs;

	Node()
		: inputs(), outputs()
		{}

	link_t * find_link_to(Node * to)
		{
		for (link_t * link : outputs)
			if (link->to == to)
				return link;
		
		return 0;
		}

	const link_t * find_link_to(const Node * to) const
		{
		for (const link_t * link : outputs)
			if (link->to == to)
				return link;
		
		return 0;
		}
	};


template<template<class> class NODE, template<class> class LINK>
struct Graph 
	{
	typedef Graph<NODE, EDGE> Self;
	typedef NODE<Self> node_t;
	typedef LINK<Self> edge_t;
	};
	

#endif	// GENERICGRAPH_H

