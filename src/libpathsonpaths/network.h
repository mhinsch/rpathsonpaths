#ifndef NETWORK_H
#define NETWORK_H

#include <vector>

using std::size_t;

/** Abstract base class for Network. */
struct AbstractNetwork
	{
	virtual void add_link(size_t from, size_t to, double rate) = 0;
	virtual void set_source(size_t s, double p) = 0;
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

	const Network & operator=(Network tmp)
		{
		swap(tmp.nodes, this->nodes);
		swap(tmp.links, this->links);
		}

	/** Add an edge. Source and target nodes have to be specified as indices.*/
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

	/** Make node @a s an external source with proportion of infected set to @a prop_infd. */
	void set_source(size_t s, double prop_infd)
		{
		assert(nodes.size() > s && nodes[s] != 0);
		
		nodes[s]->rate_in = 1.0;
		nodes[s]->rate_in_infd = prop_infd;
		nodes[s]->d_rate_in_infd = 0;
		}

	size_t find_link(L * l) const
		{
		for (size_t i=0; i<links.size(); i++)
			if (l == links[i])
				return i;

		return links.size();
		}

	size_t find_node_id(const N * n) const
		{
		for (size_t i=0; i<nodes.size(); i++)
			if (n == nodes[i])
				return i;

		return nodes.size();
		}


	~Network()
		{
		for (N * n : nodes)
			delete n;
		for (L * l : links)
			delete l;
		}

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
				assert(li < links.size());
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
				assert(li < links.size());
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
