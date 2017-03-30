#ifndef DRIFTAPPROX_H
#define DRIFTAPPROX_H

template<class NODE, class CONT>
struct DriftNode : public NODE
	{
	typedef CONT freq_t;

	freq_t frequencies;
	};

template<class NODE, class DRIFT_FUNC>
void annotate_frequencies(NODE * node, DRIFT_FUNC & drift)
	{
	assert(node->frequencies.empty());

	for (auto * link : node->inputs)
		{
		const NODE::freq_t & freq_in = link->from->frequencies;
		if (freq_in.empty())
			annotate_frequencies(link->from, drift_func);

		if (node->frequencies.empty())
			node->frequencies.resize(freq_in.size());

		const double prop = link->rate_infd / node->rate_in_infd;

		auto it_n = node->frequencies.begin();
		for (auto it_i=freq_in.begin(); it_i!=freq_in.end(); it_i++,it_n++)
			{
			}
		}
	}


#endif	// DRIFTAPPROX_H
