#ifndef DRIFTAPPROX_H
#define DRIFTAPPROX_H

#include <cassert>

template<class CONT, class NODE>
struct DriftNode : public NODE
	{
	typedef CONT freq_t;

	using NODE::NODE;

	freq_t frequencies;
	};

template<class NODE, class DRIFT_FUNC>
void annotate_frequencies(NODE * node, DRIFT_FUNC & drift)
	{
	assert(node->frequencies.empty());

	for (auto * link : node->inputs)
		{
		const typename NODE::freq_t & freq_in = link->from->frequencies;
		if (freq_in.empty())
			annotate_frequencies(link->from, drift);

		if (node->frequencies.empty())
			node->frequencies.resize(freq_in.size());

		const double prop = link->rate_infd / node->rate_in_infd;

		auto new_freq_in = freq_in;
		drift(new_freq_in);

		auto it_n = node->frequencies.begin();
		for (auto it_i=new_freq_in.begin(); it_i!=new_freq_in.end(); it_i++,it_n++)
			*it_n += *it_i * prop;
		}
	}

template<class ITER, class DRIFT_FUNC>
void annotate_frequencies(const ITER & beg, const ITER & end, DRIFT_FUNC & drift)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_frequencies(*i, drift);
	}

#endif	// DRIFTAPPROX_H
