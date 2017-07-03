#ifndef DRIFTAPPROX_H
#define DRIFTAPPROX_H

#include <numeric>


template<class NODE, class DRIFT_FUNC>
void annotate_frequencies(NODE * node, DRIFT_FUNC & drift)
	{
	if (node->done)
		return;

	if (node->is_root())
		{
		node->done = true;
		return;
		}

	// do parents in any case
	for (auto * link : node->inputs)
		annotate_frequencies(link->from, drift);

	// non-empty list indicates freqs have been pre-set 
	if (!node->frequencies.empty())
		{
		node->done = true;
		return;
		}

	// amount of incoming infected material
	const double prop_in_infd = node->rate_in_infd - node->d_rate_in_infd;

	// this is not very elegant, but I can't think of a better way to do it
	// without creating lots of little vectors all the time
	static typename NODE::freq_t res;

	for (auto * link : node->inputs)
		{
		// proportion of incoming infected material that comes from this link
		const double prop = prop_in_infd <= 0 ? 0 : link->rate_infd / prop_in_infd;

		if (prop <= 0)
			continue;

		const typename NODE::freq_t & freq_in = link->from->frequencies;
		res.resize(freq_in.size());

		drift(freq_in, res);
		// NOTE: any recursive calls will invalidate res!

		// has to go here, otherwise we won't know size
		if (node->frequencies.empty())
			node->frequencies.resize(freq_in.size(), 0);

		// doesn't look like it, but this is safe since res can never be bigger than
		// node->frequencies (unless users supply differently sized allele freqs but
		// then they are on their own)
		auto f_iter = node->frequencies.begin();
		for (const auto r : res)
			 *f_iter++ += r * prop;
		}

	node->done = true;
	}


template<class NODE, class DRIFT_FUNC>
void annotate_frequencies_push(NODE * node, DRIFT_FUNC & drift)
	{
	if (node->done)
		return;

	// do parents in any case
	for (auto l : node->inputs)
		annotate_frequencies_push(l->from, drift);

	// we want even empty nodes to have a set of frequencies
	// so let's do that here
	for (auto l : node->outputs)
		l->to->frequencies.resize(node->frequencies.size(), 0.0);
	
	// we are pushing, so ignore leaves
	if (node->is_leaf() || node->rate_in <= 0 || node->rate_in_infd <= 0)
		{
		node->done = true;
		return;
		}

	// this branch of the graph is dead
	if (node->frequencies.empty())
		{
		node->done = true;
		return;
		}

	// this is not very elegant, but I can't think of a better way to do it
	// without creating lots of little vectors all the time
	static typename NODE::freq_t res;
	res.resize(node->frequencies.size());

	for (auto link : node->outputs)
		{
		auto to = link->to;
		// pre-transmission infected in target node
		const double to_n_infd = to->rate_in_infd - to->d_rate_in_infd;

		if (to_n_infd <= 0) continue;

		// proportion of those coming through this link
		const double p_to = link->rate_infd / to_n_infd;

		// link rate might be zero
		if (p_to <= 0) continue;

		drift(node->frequencies, res);
		// NOTE: any recursive calls will invalidate res!

		// doesn't look like it, but this is safe since res can never be bigger than
		// node->frequencies (unless users supply differently sized allele freqs but
		// then they are on their own)
		auto f_iter = to->frequencies.begin();
		for (const auto r : res)
			 *f_iter++ += r * p_to;
		}

	node->done = true;
	}


template<class ITER, class DRIFT_FUNC>
void annotate_frequencies(const ITER & beg, const ITER & end, DRIFT_FUNC & drift)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_frequencies_push(*i, drift);

	for (ITER i=beg; i!=end; i++)
		(*i)->done = false;
	}


#endif	// DRIFTAPPROX_H
