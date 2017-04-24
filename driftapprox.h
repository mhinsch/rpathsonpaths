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
		{
		assert(!node->frequencies.empty());
		return;
		}

	assert(node->frequencies.empty());

	const double prop_in_infd = node->rate_in_infd - node->d_rate_in_infd;

	for (auto * link : node->inputs)
		{
		typename NODE::freq_t & freq_in = link->from->frequencies;
		if (freq_in.empty())
			annotate_frequencies(link->from, drift);

		if (node->frequencies.empty())
			node->frequencies.resize(freq_in.size(), 0);

		const double prop = link->rate_infd / prop_in_infd;

		drift(freq_in);

		auto it_i = drift.begin();
		for (auto & freq : node->frequencies)
			{
			freq += *it_i * prop;
			it_i++;
			}
		}
	}

template<class ITER, class DRIFT_FUNC>
void annotate_frequencies(const ITER & beg, const ITER & end, DRIFT_FUNC & drift)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_frequencies(*i, drift);
	}

#endif	// DRIFTAPPROX_H
