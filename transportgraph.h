#ifndef TRANSPORTGRAPH_H
#define TRANSPORTGRAPH_H

#include <vector>

using namespace std;

template<class LINK>
struct TranspLink : public LINK
	{
	double rate;
	double rate_infd;

	TranspLink(LINK::node_t * f = 0, LINK::node_t * t = 0, double r = 0.0, double ir = -1)
		: LINK(f, t), rate(r), rate_infd(ir)
		{}

	bool is_source() const
		{
		return from == 0;
		}
	};


template<class NODE>
struct TranspNode : public NODE
	{
	double rate_in, rate_in_infd, d_rate_in_infd;
	double rate_out_infd;

	TranspNode()
		: NODE(), rate_in(0), rate_in_infd(0), rate_out_infd(0)
		{}

	// probability an infected unit coming from this node was newly infected
	double prob_newly_infected() const
		{
		// delta inf / inf
		return  d_rate_in_infd / rate_in_infd;
		}
	};

/* current implicit assumption:
 * - sources are dangling input links with rate_infd already set
 * - links: rate_infd < 0 <=> haven't been processed yet
 * - nodes: rate_in <= 0 <=> haven't been processed yet
 */

template<class NODE>
void annotate_rates(NODE * node, double transm_rate)
	{
	// this one is done
	if (node->rate_in > 0)
		return;

	// *** input

	for (NODE::link_t * link : node->inputs)
		{
		if (link->rate_infd < 0)
			annotate_rates(link->from);

		node->rate_in += link->rate;
		node->rate_in_infd += link->rate_infd;
		}

	// *** infection
	
	// proportion of input becomes infected
	node->d_rate_in_infd = transm_rate * (node->rate_in - node->rate_in_infd);
	node->rate_in_infd += node->d_rate_infd;

	// proportion of infected units
	double prop_infd = node->rate_in_infd / node->rate_in;	
	
	// *** output

	for (NODE::link_t * link : node.outputs)
		{
		// if this has been done there's something wrong
		assert(link->rate_infd < 0);
		
		// all outputs have the same proportion of infected units
		// NOTE could be changed down the line
		link->rate_infd = link->rate * prop_infd;

		node->rate_out_infd += link->rate_infd;
		}
	}

template<class ITER>
void annotate_rates(const ITER & beg, const ITER & end, double transm_rate)
	{
	for (ITER i=beg, i!=end; i++)
		annotate_rates(*i, transm_rate);
	}

// probability of infected material from n_from to end up in n_to
template<class NODE>
double prob(NODE * n_from, NODE * n_to) const
	{
	const NODE::link_t * link = n_from->find_link_to(n_to);
	if (link)
		return link->rate_infd / node->rate_out_infd;
	
	assert(false);
	return 0;
	}
};


#endif	// TRANSPORTGRAPH_H
