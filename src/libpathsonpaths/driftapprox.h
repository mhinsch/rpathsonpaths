#ifndef DRIFTAPPROX_H
#define DRIFTAPPROX_H

#include <numeric>

template<class CONT>
struct FreqNode
	{
	typedef CONT freq_t;
	typedef typename freq_t::value_type value_t;

	freq_t frequencies;

	bool valid(value_t epsilon) const
		{
		const auto sum = std::accumulate(frequencies.begin(), frequencies.end(), value_t(0));
		return sum + epsilon > 1 && sum - epsilon < 1;
		}

	value_t normalize(value_t norm=1.0)
		{
		const auto sum = std::accumulate(frequencies.begin(), frequencies.end(), value_t(0));
		if (sum > 0 && sum != norm)
			{
			const auto factor = norm / sum;
			for (auto & f : frequencies)
				f *= factor;
			}

		return sum;
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

	// amount of incoming infected material
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

		// proportion of incoming infected material that comes from this link
		const double prop = prop_in_infd <= 0 ? 0 : link->rate_infd / prop_in_infd;

		// doesn't look like it, but this is safe since res can never be bigger than
		// node->frequencies (unless users supply differently sized allele freqs but
		// then they are on their own)
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
