#ifndef GENEFREQGRAPH_H
#define GENEFREQGRAPH_H

#include <numeric>

/** A node class that keeps a list of allele frequencies. */
template<class CONT>
struct FreqNode
	{
	typedef CONT freq_t;
	typedef typename freq_t::value_type value_t;

	freq_t frequencies; //!< allele frequencies

	/** Check whether frequencies some up to 1 (+- epsilon). */
	bool valid(value_t epsilon) const
		{
		const auto sum = std::accumulate(frequencies.begin(), frequencies.end(), value_t(0));
		return sum + epsilon > 1 && sum - epsilon < 1;
		}

	/** Rescale frequencies to a given sum. */
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


#endif	// GENEFREQGRAPH_H
