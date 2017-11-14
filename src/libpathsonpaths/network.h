#ifndef NETWORK_H
#define NETWORK_H

#include <vector>

#include "util.h"

using std::size_t;

/** Abstract base class for Network. */
struct AbstractNetwork
	{
	/** Add a link to the network. 
	 * @param from Start node.
	 * @param to End node.
	 * @param rate Transport rate. */
	virtual void add_link(size_t from, size_t to, double rate) = 0;
	/** Mark a node as source. 
	 * @param s Source node.
	 * @param p Proportion of infected material.
	 * @param i Overall input at this node. */
	virtual void set_source(size_t s, double p, double i=1.0) = 0;
	/** Plain destructor. */
	virtual ~AbstractNetwork()
		{
		}
	};


/** Network template class. Network keeps track of nodes and edges and provides methods
 * to construct a network by adding links. Note that Node/Link objects are assumed to be 
 * owned by Network.
 * @param N node type.
 * @param L link type. */
template<class N, class L>
struct Network : public AbstractNetwork
	{
	std::vector<N *> nodes;		//!< all nodes in the network
	std::vector<L *> links;		//!< all edges in the network

	Network() = default;

	Network(const Network & other)
		{
		other.clone_into(*this);
		}

	const Network & operator=(Network && tmp)
		{
		swap(tmp.nodes, this->nodes);
		swap(tmp.links, this->links);

		return *this;
		}

	/** Add an edge. Source and target nodes have to be specified as indices.
	 * @param from, to indices of nodes being linked
	 * @param transfer rate of material (in amount per unit time) */
	void add_link(size_t from, size_t to, double rate)
		{
		if (nodes.size() <= std::max(from, to))
			nodes.resize(std::max(from, to)+1, 0);
		if (nodes[from] == 0)
			nodes[from] = new N;
		if (nodes[to] == 0)
			nodes[to] = new N;

		links.push_back(new L(nodes[from], nodes[to], rate));
		nodes[from]->add_output(links.back());
		nodes[to]->add_input(links.back());
		}

	void set_source(size_t s, double p, double i) {}

	/** Find index of link @a l. */
	size_t find_link(L * l) const
		{
		for (size_t i=0; i<links.size(); i++)
			if (l == links[i])
				return i;

		return links.size();
		}

	/** Find index of node @a n. */
	size_t find_node_id(const N * n) const
		{
		for (size_t i=0; i<nodes.size(); i++)
			if (n == nodes[i])
				return i;

		return nodes.size();
		}

	/** Reset done status to false for all nodes. */
	void reset_done()
		{
		for (auto n : nodes)
			n->done = false;
		}

	/** Destructor. Deletes all nodes and links. */
	~Network()
		{
		for (N * n : nodes)
			delete n;
		for (L * l : links)
			delete l;
		}

	/** Deep copy the content of this network into another one. Any previous content in the 
	 * other network will be lost. */
	void clone_into(Network & nn) const
		{
		nn.nodes.resize(nodes.size(), 0);
		nn.links.resize(links.size(), 0);

		// copy all links, pointer will be readjusted later
		for (size_t i=0; i<links.size(); i++)
			nn.links[i] = new L(*links[i]);

		for (size_t i=0; i<nodes.size(); i++)
			{
			// copy node (pointers and all)
			N * n = new N(*nodes[i]);

			// re-point input pointers
			for (auto & l : n->inputs)
				{
				// this works because the node still has the old pointers
				size_t li = find_link(l);
				myassert(li < links.size());
				// use new link object
				l = nn.links[li];
				// point it to this node
				l->to = n;
				}

			// re-point output pointers
			for (auto & l : n->outputs)
				{
				// this works because the node still has the old pointers
				size_t li = find_link(l);
				myassert(li < links.size());
				// use new link object
				l = nn.links[li];
				// point it to this node
				l->from = n;
				}
			
			nn.nodes[i] = n;
			}
		}
	};


#endif	// NETWORK_H
