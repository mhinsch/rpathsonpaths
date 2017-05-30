#ifndef TRANSPORTGRAPH_H
#define TRANSPORTGRAPH_H

#include "util.h"

/** A link type that keeps track of transfer rates.  */
struct TranspLink 
	{
	double rate;		//!< (absolute) transfer rate of material.
	double rate_infd;	//!< (absolute) transfer rate of infected material.

	TranspLink(double a_rate = 0.0, double a_rate_infd = -1)
		: rate(a_rate), rate_infd(a_rate_infd)
		{}
	};


/** A node type that keeps track of input and output rates.*/
struct TranspNode 
	{
	double rate_in;			//!< overall input rate.
	double rate_in_infd;	//!< overall rate of input of infected material (@a after transmission).
	double d_rate_in_infd;	//!< 
	double rate_out_infd;

	TranspNode()
		: rate_in(-1), rate_in_infd(0), d_rate_in_infd(0), rate_out_infd(-1)
		{}

	void reset_rates()
		{
		rate_in = 0;
		rate_in_infd = 0;
		d_rate_in_infd = 0;
		rate_out_infd = 0;
		}

	double prop_infected() const
		{
		return rate_in <= 0 ?  
			0 : rate_in_infd / rate_in;	
		}

	/** Probability an infected unit coming from this node was newly infected. */
	double prob_newly_infected() const
		{
		// delta inf / inf
		return rate_in_infd <= 0 ? 0 : d_rate_in_infd / rate_in_infd;
		}
	};


/** Adjust output rates so that sum(output) = sum(input) * (1-decay). */
template<class NODE>
void preserve_mass(NODE * node, double decay)
	{
	if (node->done || node->is_leaf())
		return;

	for (auto i : node->inputs)
		preserve_mass(i->from, decay);

	double inp = 0.0;

	// root nodes use preset value
	// we could allow for more flexibility (e.g. use preset if inp == 0)
	// but this is the most consistent option
	if (node->is_root())
		inp = node->rate_in;
	else	
		for (auto l : node->inputs)
			inp += l->rate;

	double outp = 0.0;
	for (auto l : node->outputs)
		outp += l->rate;

	myassert(outp > 0);

	const double f = (inp * (1.0 - decay)) / outp;

	for (auto l : node->outputs)
		l->rate *= f;

	node->done = true;
	}


template<class ITER>
void preserve_mass(const ITER & beg, const ITER & end, double decay)
	{
	for (ITER i=beg; i!=end; i++)
		(*i)->done = false;

	for (ITER i=beg; i!=end; i++)
		preserve_mass(*i, decay);
	}


/** Calculate overall rate of infected input and proportion of infected material
 * in NODE node (after transmission) and in its output. 
 *
 * @note This function will call itself 
 * recursively for unprocessed input nodes.
 *
 * @note Current implicit assumption:
 * - sources are input-less nodes with rates pre-set
 * - links: rate_infd < 0 <=> haven't been processed yet
 * - nodes: rate_in <= 0 <=> haven't been processed yet
 * 
 * @tparam NODE node type.
 * @param node node to process.
 * @param transm_rate rate of infection within nodes
 */
template<class NODE>
void annotate_rates(NODE * node, double transm_rate)
	{
	// not processed yet
	// NOTE sources will have this set but *not* the output rate!
	if (node->rate_in < 0)
		{
		// *** input
		
		node->rate_in = 0;

		// does nothing for roots
		for (typename NODE::link_t * link : node->inputs)
			{
			// new links set that to -1
			if (link->rate_infd < 0)
				annotate_rates(link->from, transm_rate);

			node->rate_in += link->rate;
			node->rate_in_infd += link->rate_infd;
			}

		// *** infection
		
		// proportion of input becomes infected
		node->d_rate_in_infd = transm_rate * (node->rate_in - node->rate_in_infd);
		node->rate_in_infd += node->d_rate_in_infd;
		}

	if (node->rate_out_infd < 0)
		{
		// *** output
		
		node->rate_out_infd = 0;
		
		// proportion of infected units
		const double prop_infd = node->prop_infected();

		// does nothing for leaves
		for (typename NODE::link_t * link : node->outputs)
			{
			// if this has been done there's something wrong
			myassert(link->rate_infd < 0);
			
			// all outputs have the same proportion of infected units
			// NOTE could be changed down the line
			link->rate_infd = link->rate * prop_infd;

			node->rate_out_infd += link->rate_infd;
			}
		}
	}

/** Annotate rates for a collection of nodes and all their inputs (recursively). 
 *
 * @tparam ITER iterator over nodes.
 * @param beg, end range of nodes to be processed.
 * @param transm_rate infection rate within nodes.
 */
template<class ITER>
void annotate_rates(const ITER & beg, const ITER & end, double transm_rate)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_rates(*i, transm_rate);
	}

/** Probability of infected material from node @a n_from to end up in node @a n_to. 
 *
 * @pre Assumes that there is a link from @a n_from to @a n_to.
 * 
 * @tparam NODE node type.
 */
template<class NODE>
double prob(NODE * n_from, NODE * n_to) 
	{
	const typename NODE::link_t * link = n_from->find_link_to(n_to);
	if (link)
		return n_from->rate_out_inf <= 0 ? 0 : link->rate_infd / n_from->rate_out_infd;
	
	myassert(false);
	return 0;
	}



#endif	// TRANSPORTGRAPH_H
