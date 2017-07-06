#ifndef GENERICGRAPH_H
#define GENERICGRAPH_H

/** Generic Link type.
 * @tparam GRAPH helper class that provides the node type.
 */
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

/** Generic Node type. Provides simple helper functions related to network topology.
 * @tparam GRAPH helper class that provides the link type.
 * @tparam CONT container type to be used to store links.
 */
template<class GRAPH, template<class> class CONT>
struct Node
	{
	typedef GRAPH graph_t;
	typedef typename graph_t::link_t link_t;

	typedef CONT<link_t *> cont_t;

	cont_t inputs;
	cont_t outputs;

	bool done;

	Node()
		: inputs(), outputs(), done(false)
		{}

	void add_input(link_t * inp)
		{
		inputs.push_back(inp);
		}
	void add_output(link_t * outp)
		{
		outputs.push_back(outp);
		}

	/** Find a link in container @a c that connects to Node @a to.
	 * @tparam FORWARD whether the link is a forward or a backward link.
	 */
	template<bool FORWARD = true>
	static link_t * find_link(Node * to, cont_t & c)
		{
		for (link_t * link : c)
			if ((FORWARD ? link->to : link->from) == to)
				return link;
		
		return 0;
		}

	template<bool FORWARD = true>
	static const link_t * find_link(const Node * to, const cont_t & c)
		{
		for (const link_t * link : c)
			if ((FORWARD ? link->to : link->from) == to)
				return link;
		
		return 0;
		}

	link_t * find_link_to(Node * to)
		{
		return find_link<true>(to, outputs);
		}
	const link_t * find_link_to(const Node * to) const
		{
		return find_link<true>(to, outputs);
		}

	link_t * find_link_from(Node * from)
		{
		return find_link<false>(from, inputs);
		}
	const link_t * find_link_from(const Node * from) const
		{
		return find_link<false>(from, inputs);
		}

	/** Check for network consistency. In particular this checks whether all
	 * input nodes have the current node as an output and whether all output
	 * nodes have the current node as an input.
	 */
	bool consistent() const 
		{ 
		for (const link_t * link : inputs) 
			{ 
			const link_t * l = link->from->find_link_to(this); 
			if (!l || l!=link || l->to != this) 
				return false; 
			}

		for (const link_t * link : outputs)
			{
			const link_t * l = link->to->find_link_from(this);
			if (!l || l!=link || l->from != this)
				return false;
			}

		return true;
		}

	/** A leaf node is a node with no outputs. */
	bool is_leaf() const
		{
		return outputs.size() == 0;
		}
	/** A root node is a node with no inputs. */
	bool is_root() const
		{
		return inputs.size() == 0;
		}
	};

/** Helper type binding a node and a link type together. */
template<template<class> class NODE, template<class> class LINK>
struct Graph 
	{
	typedef Graph<NODE, LINK> Self;
	typedef NODE<Self> node_t;
	typedef LINK<Self> link_t;
	};
	

template<class NODE>
void reset_downstream(NODE & node, bool to=false)
	{
	node.done = to;

	for (auto l : node.outputs)
		reset_downstream(*l->to, to);
	}


template<class NODE>
void reset_upstream(NODE & node, bool to=false)
	{
	node.done = to;

	for (auto l : node.inputs)
		reset_upstream(*l->from, to);
	}

template<bool PREORDER=true, class NODE, class FUNC>
void apply_downstream(NODE & node, FUNC func)
	{
	if (node.done)
		return;

	node.done = true;

	if (PREORDER)
		func(node);

	for (auto l : node.outputs)
		apply_downstream(*l->to, func);

	if (!PREORDER)
		func(node);
	}


template<bool PREORDER=true, class NODE, class FUNC>
void apply_upstream(NODE & node, FUNC func)
	{
	if (node.done)
		return;

	node.done = true;

	if (PREORDER)
		func(node);

	for (auto l : node.inputs)
		apply_upstream(*l->from, func);

	if (!PREORDER)
		func(node);
	}

#endif	// GENERICGRAPH_H

