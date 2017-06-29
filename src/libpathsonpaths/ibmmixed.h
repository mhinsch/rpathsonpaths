#ifndef IBMMIXED_H
#define IBMMIXED_H

#include <vector>

#include "util.h"

template<class NODE, class RNG>
void annotate_rates_ibmm(NODE * node, double transm_rate, RNG & rng)
	{
	if (node->done)
		return;

	// not processed yet
	// NOTE sources will have this set but *not* the output rate!
	if (node->rate_in < 0)
		{
		node->rate_in = 0;

// *** collect input

		// this needs to be done first to make sure we get the whole network set up
		// properly
		for (auto link : node->inputs)
			{
			annotate_rates_ibmm(link->from, transm_rate, rng);

			node->rate_in += link->rate;
			node->rate_in_infd += link->rate_infd;
			}

// *** transmission
//     NOTE that this doesn't happen for sources

		const double inp = node->rate_in;
		const double uninfd = inp - node->rate_in_infd;
		// coin flip for each uninfected whether it becomes infected
		const int newly_infd = uninfd > 0 ? rng.binom(transm_rate, int(uninfd)) : 0;

		node->rate_in_infd += newly_infd;
		node->d_rate_in_infd = newly_infd;

		myassert(inp >= 0 && uninfd >= 0 && newly_infd >=0);
		}

// *** output

	double outp = 0.0;
	for (const auto & l : node->outputs)
		{
		outp += l->rate;
		// do that here so that we can just quit if there's no output
		l->rate_infd = 0;
		}

	node->rate_out_infd = 0;

	// no output, done
	if (outp <= 0)
		{
		node->done = true;
		return;
		}
	
// *** generate output
//
// We don't replace units that have been selected for output => we have to 
// use a hypergeometric distribution 	

	// we have to keep track of how many infected there are still left
	// all_infd == sum(node->frequencies)
	double all_infd = node->rate_in_infd;
	double all_non_infd = node->rate_in - node->rate_in_infd;
	myassert(all_non_infd >= 0);

	//std::cout << "****\n";
	for (auto l : node->outputs)
		{
		// how many go into this link (infd + non-infd)
		// needs to be an int, otherwise we'll get into trouble b/c rounding errors
		int pick = int(l->rate);

		myassert(pick <= all_infd + all_non_infd);

		const double a = rng.hypergeom(int(all_infd), int(all_non_infd), int(pick));
		l->rate_infd = a;
		all_infd -= a;
		all_non_infd -= (pick - a);
		}

	// adjust output rates in the node
	for (const auto & l : node->outputs)
		node->rate_out_infd += l->rate_infd;

	node->done = true;
	}


template<class ITER, class RNG>
void annotate_rates_ibmm(const ITER & beg, const ITER & end, double transm_rate, RNG & rng)
	{
	for (ITER i=beg; i!=end; i++)
		(*i)->done = false;

	for (ITER i=beg; i!=end; i++)
		annotate_rates_ibmm(*i, transm_rate, rng);
	}


template<class NODE, class RNG>
void annotate_frequencies_ibmm(NODE * node, RNG & rng)
	{
	if (node->done)
		return;

	// this needs to be done first to make sure we get the whole network set up
	// properly
	for (auto link : node->inputs)
		annotate_frequencies_ibmm(link->from, rng);

	// we want even empty nodes to have a set of frequencies
	// so let's do that here
	for (auto l : node->outputs)
		l->to->frequencies.resize(node->frequencies.size(), 0.0);

	// we are pushing, so ignore leaves
	if (node->is_leaf() || node->rate_in <= 0)
		{
		node->normalize();
		node->done = true;
		return;
		}

	// no input set on this branch
	if (node->frequencies.size() == 0)
		{
		node->done = true;
		return;
		}

	if (node->is_root())
		{
		// just in case, frequencies might have been initialized to a different scale than
		// input rate; we just assume input rate has priority
		node->normalize(node->rate_in_infd);
		// now we move everything to int, otherwise bad things will happen later on
		int sum = 0;
		for (auto & f : node->frequencies)
			sum += (f = int(f));
		node->rate_in_infd = sum;
		}

	double outp = 0.0;
	for (const auto & l : node->outputs)
		outp += l->rate;

	if (outp <= 0)
		{
		node->normalize();
		node->done = true;
		return;
		}

	// pre-transmission infected
	const double infd = node->rate_in_infd - node->d_rate_in_infd;

	if (infd <= 0)
		{
		node->done = true;
		return;
		}

	const double inp = node->rate_in;
	const int newly_infd = int(node->d_rate_in_infd);

	myassert(inp > 0 && newly_infd >=0);
	
// *** transmission 
//
// Units that transmitted infection remain in the pool => polynomial distribution
// (or binomial + conditional method).

	// if there are newly infected units, distribute them according to prop. of 
	// alleles; has to check for infected as well (due to weird infection model)
	// multinomial assumes that all infections happen simultaneously
	if (newly_infd > 0)
		{
		int n = newly_infd;
		double rem = 1.0;

		for (size_t i=0; i<node->frequencies.size()-1 && n>0 && rem>0; i++)
			{
			const double p = node->frequencies[i] / infd;
			// due to numeric effects it can happen that p>rem (slightly)
			// if this is the last positive frequency
			const int add = rng.binom(std::min(1.0, p/rem), n);

			myassert(add >= 0);

			node->frequencies[i] += add;

			n -= add;
			rem -= p;
			}

		myassert(n>=0 && rem>-0.0001); 
		// R binom does weird stuff when rem and p are very close so we
		// skip the last step and just assign n directly
		node->frequencies.back() += n;
		}

// *** generate output
//
// We don't replace units that have been selected for output => we have to 
// use a multi-variate hypergeometric distribution (we actually use a regular 
// hypergeometric and the conditional method).

	// we have to keep track of how many infected there are still left
	// left_all == sum(node->frequencies)
	double left_all = node->rate_in_infd;
	auto left_by_gene = node->frequencies;

	for (auto l : node->outputs)
		{
		// how many infd go into this link (uninfd have been done by annotate_rates)
		// needs to be an int, otherwise we'll get into trouble b/c rounding errors
		int pick = int(l->rate_infd);
		// how many to pick from for the remaining alleles
		double all_infd = left_all;

		// first n-1 alleles
		for (size_t i=0; i<left_by_gene.size()-1; i++)
			{
			// split leftover into two groups
			all_infd -= left_by_gene[i];
			
			myassert(all_infd >= 0);
			myassert(pick <= all_infd+left_by_gene[i]);

			// picked units are either allele i or allele >i
			const int add = rng.hypergeom(int(left_by_gene[i]), int(all_infd), int(pick));
			
			myassert(add >= 0);
			// adjust number to pick from for the next output
			left_by_gene[i] -= add;
			myassert(left_by_gene[i] >= 0);
			// keep this synchronized with left_by_gene
			left_all -= add;
			myassert(left_all >= 0);
			// we pick less next time
			pick -= add;
			myassert(pick >= 0);
			// add the stuff to the link and its end node
			l->to->frequencies[i] += add;
			}

		// last allele gets the leftovers
		l->to->frequencies.back() += pick;
		// keep track
		left_all -= pick;
		left_by_gene.back() -= pick;
		}

	node->normalize();
	node->done = true;
	}


template<class ITER, class BINOM_FUNC>
void annotate_frequencies_ibmm(const ITER & beg, const ITER & end, BINOM_FUNC & binom)
	{
	for (ITER i=beg; i!=end; i++)
		(*i)->done = false;

	for (ITER i=beg; i!=end; i++)
		annotate_frequencies_ibmm(*i, binom);
	}


#endif	// IBMMIXED_H
