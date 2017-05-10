#ifndef DRIFTAPPROX_H
#define DRIFTAPPROX_H

#include <cassert>

template<class CONT, class NODE>
struct DriftNode : public NODE
	{
	typedef CONT freq_t;

	using NODE::NODE;

	freq_t frequencies;

	bool valid(typename freq_t::value_type epsilon) const
		{
		typename freq_t::value_type sum(0);
		for (auto f : frequencies)
			sum += f;

		return sum + epsilon > 1 && sum - epsilon < 1;
		}
	};

template<class NODE, class DRIFT_FUNC>
void annotate_frequencies(NODE * node, DRIFT_FUNC & drift)
	{
	if (node->is_root())
		return;

	// non-empty list indicates already processed
	if (!node->frequencies.empty())
		return;

	const double prop_in_infd = node->rate_in_infd - node->d_rate_in_infd;

	// this is not very elegant, but I can't think of a better way to do it
	// without creating lots of little vectors all the time
	static typename NODE::freq_t res;

	for (auto * link : node->inputs)
		{
		if (link->from->frequencies.empty())
			annotate_frequencies(link->from, drift);

		const typename NODE::freq_t & freq_in = link->from->frequencies;
		res.resize(freq_in.size());

		drift(freq_in, res);
		// NOTE: any recursive calls will invalidate res!

		// has to go here, otherwise we won't know size
		if (node->frequencies.empty())
			node->frequencies.resize(freq_in.size(), 0);

		const double prop = prop_in_infd <= 0 ? 0 : link->rate_infd / prop_in_infd;

		auto f_iter = node->frequencies.begin();
		for (const auto r : res)
			 *f_iter++ += r * prop;
		}
	}

template<class ITER, class DRIFT_FUNC>
void annotate_frequencies(const ITER & beg, const ITER & end, DRIFT_FUNC & drift)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_frequencies(*i, drift);
	}

#endif	// DRIFTAPPROX_H
