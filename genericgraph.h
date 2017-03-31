#ifndef GENERICGRAPH_H
#define GENERICGRAPH_H

template<class GRAPH>
struct Link
	{
	typedef GRAPH graph_t;
	typedef typename graph_t::node_t node_t;

	node_t * from, * to;

	Link(node_t * f = 0, node_t * t = 0)
		: from(f), to(t)
		{}
	};

template<class GRAPH, template<class> class CONT>
struct Node
	{
	typedef GRAPH graph_t;
	typedef typename graph_t::link_t link_t;

	typedef CONT<link_t *> cont_t;

	cont_t inputs;
	cont_t outputs;

	Node()
		: inputs(), outputs()
		{}

	void add_input(link_t * inp)
		{
		inputs.push_back(inp);
		}
	void add_output(link_t * outp)
		{
		outputs.push_back(outp);
		}

	static link_t * find_link(Node * to, cont_t & c)
		{
		for (link_t * link : c)
			if (link->to == to)
				return link;
		
		return 0;
		}

	static const link_t * find_link(const Node * to, const cont_t & c)
		{
		for (const link_t * link : c)
			if (link->to == to)
				return link;
		
		return 0;
		}

	link_t * find_link_to(Node * to)
		{
		return find_link(to, outputs);
		}
	const link_t * find_link_to(const Node * to) const
		{
		return find_link(to, outputs);
		}

	link_t * find_link_from(Node * from)
		{
		return find_link(from, inputs);
		}
	const link_t * find_link_from(const Node * from) const
		{
		return find_link(from, inputs);
		}

	bool is_leaf() const
		{
		return outputs.size() == 0;
		}
	bool is_root() const
		{
		return inputs.size() == 0;
		}
	};


template<template<class> class NODE, template<class> class LINK>
struct Graph 
	{
	typedef Graph<NODE, LINK> Self;
	typedef NODE<Self> node_t;
	typedef LINK<Self> link_t;
	};
	

#endif	// GENERICGRAPH_H

