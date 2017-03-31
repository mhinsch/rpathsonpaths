#ifndef NETWORK_H
#define NETWORK_H

#include <vector>

using std::size_t;

struct AbstractNetwork
	{
	virtual void add_link(size_t from, size_t to, double rate) = 0;
	virtual void set_source(size_t s, double p) = 0;
	};


template<class N, class L>
struct Network : public AbstractNetwork
	{
	std::vector<N *> nodes;
	std::vector<L *> links;

	void add_link(size_t from, size_t to, double rate)
		{
		if (nodes.size() <= std::max(from, to))
			nodes.resize(std::max(from, to));
		if (nodes[from] == 0)
			nodes[from] = new N;
		if (nodes[to] == 0)
			nodes[to] = new N;

		links.push_back(new L(nodes[from], nodes[to], rate));
		nodes[from]->add_output(links.back());
		nodes[to]->add_input(links.back());
		}

	void set_source(size_t s, double prop_infd)
		{
		assert(nodes.size() > s && nodes[s] != 0);
		
		nodes[s]->rate_in = 1.0;
		nodes[s]->rate_in_infd = prop_infd;
		nodes[s]->d_rate_in_infd = 0;
		}

	~Network()
		{
		for (N * n : nodes)
			delete n;
		for (L * l : links)
			delete l;
		}
	};


#endif	// NETWORK_H
