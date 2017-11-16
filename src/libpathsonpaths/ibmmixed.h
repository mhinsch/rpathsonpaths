#ifndef IBMMIXED_H
#define IBMMIXED_H

#include <vector>

#include "util.h"

/** Run mechanistic infection and spread simulation on node (and ancestors). */
template<class NODE, class RNG>
void annotate_rates_ibmm(NODE * node, double transm_rate, RNG & rng)
	{
	if (node->done)
		return;

// *** collect input

	if (!node->is_root())
		node->rate_in = node->rate_in_infd = 0;

	// this needs to be done first to make sure we get the whole network set up
	// properly
	for (auto link : node->inputs)
		{
		annotate_rates_ibmm(link->from, transm_rate, rng);

		node->rate_in += link->rate;
		node->rate_in_infd += link->rate_infd;
		}

	if (node->rate_in_infd <= 0)
		{
		node->done = true;
		return;
		}

	const int in_infd = node->rate_in_infd;
	//const int inp = node->rate_in;

// *** transmission
//     NOTE that this also happens for sources

	const int uninfd = node->rate_in - node->rate_in_infd;
	// coin flip for each uninfected whether it becomes infected
	// we need to check for int(infected input) otherwise we might get
	// output but 0 allele frequencies
	const int newly_infd = uninfd>0 && in_infd>0 ? 
		rng.binom(transm_rate, uninfd) : 0;

	node->rate_in_infd = in_infd + newly_infd;
	node->d_rate_in_infd = newly_infd;

	ensure(uninfd >= 0, "transport rate smaller than number of infected");
	ensure(newly_infd >=0, "negative number of new infections");

// *** output

	node->rate_out_infd = 0;

	double outp = 0.0;
	for (const auto & l : node->outputs)
		outp += l->rate;

	ensure(outp <= node->rate_in, "output can't be bigger than input");

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
	int all_infd = node->rate_in_infd;
	int all_non_infd = node->rate_in - node->rate_in_infd;
	myassert(all_non_infd >= 0);

	//std::cout << "****\n";
	for (auto l : node->outputs)
		{
		// how many go into this link (infd + non-infd)
		// needs to be an int, otherwise we'll get into trouble b/c rounding errors
		int pick = int(l->rate);

		myassert(pick <= all_infd + all_non_infd);

		const int a = rng.hypergeom(all_infd, all_non_infd, pick);
		l->rate_infd = a;
		all_infd -= a;
		all_non_infd -= (pick - a);
		}

	// adjust output rates in the node
	for (const auto & l : node->outputs)
		node->rate_out_infd += l->rate_infd;

	node->done = true;
	}


/** Run mechanistic infection and spread simulation on a range of nodes. */
template<class ITER, class RNG>
void annotate_rates_ibmm(const ITER & beg, const ITER & end, double transm_rate, RNG & rng)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_rates_ibmm(*i, transm_rate, rng);

	for (ITER i=beg; i!=end; i++)
		(*i)->done = false;
	}


/** Stochastically scale from frequency of infected individuals to absolute numbers. */
template<class NODE, class RNG>
void freq_to_popsize_ibmm(NODE * node, RNG & rng)
	{
	if (node->frequencies.empty()) 
		return;

	int n = int(node->rate_in_infd - node->d_rate_in_infd);

	if (n <= 0)
		{
		std::fill(node->frequencies.begin(), node->frequencies.end(), 0.0);
		return;
		}

	double rem = std::accumulate(node->frequencies.begin(), node->frequencies.end(), 0.0);

	ensure(rem >= 0, "negative number of infected units");

	// invalid or already scaled
	// (we have to special case 1 since that's the canonical sum of frequencies)
	if (rem <= 0 || (n>1 && rem == n))
		return;

	for (size_t i=0; i<node->frequencies.size()-1; i++)
		{
		const double p = node->frequencies[i];
		// due to numeric effects it can happen that p>rem (slightly)
		// if this is the last positive frequency
		const int add = n>0 && rem >0 ? rng.binom(std::min(1.0, p/rem), n) : 0;

		ensure(add >= 0, "internal error while scaling frequencies");

		node->frequencies[i] = add;

		n -= add;
		rem -= p;
		}

	ensure(n>=0 && rem>-0.0001, "internal error while scaling frequencies"); 
	// R binom does weird stuff when rem and p are very close so we
	// skip the last step and just assign n directly
	node->frequencies.back() = n;
	}


/** Scale frequencies to absolute numbers for a range of nodes. */
template<class ITER, class BINOM_FUNC>
void freq_to_popsize_ibmm(const ITER & beg, const ITER & end, BINOM_FUNC & binom)
	{
	// get #individuals from frequencies
	for (ITER i=beg; i!=end; i++)
		freq_to_popsize_ibmm(*i, binom);
	}

/** Run mechanistic genetics simulation. */
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
		node->done = true;
		return;
		}

	// no input set on this branch
	if (node->frequencies.empty())
		{
		node->done = true;
		return;
		}

	double outp = 0.0;
	for (const auto & l : node->outputs)
		outp += l->rate;

	// no output, done
	if (outp <= 0)
		{
		node->done = true;
		return;
		}

	ensure(outp <= node->rate_in, "output can't be bigger than input");

	// pre-transmission infected
	const double infd = node->rate_in_infd - node->d_rate_in_infd;

	myassert(std::accumulate(node->frequencies.begin(), node->frequencies.end(), 0.0) == infd);

	// nothing infected, done
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
		// tracks sum(freq[i..n])
		int infd_left = infd;

		for (size_t i=0; i<node->frequencies.size()-1; i++)
			{
			const double p = node->frequencies[i] / infd_left;
			// due to numeric effects it can happen that p>rem (slightly)
			// if this is the last positive frequency
			const int add = n>0 ? rng.binom(std::min(1.0, p), n) : 0;

			myassert(add >= 0);

			infd_left -= node->frequencies[i];
			node->frequencies[i] += add;
			n -= add;
			}

		myassert(n>=0); 
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

		// link rate might be 0
		if (pick == 0) continue;
		myassert(pick > 0);

		// how many to pick from for the remaining alleles
		double all_infd = left_all;

		// first n-1 alleles
		for (size_t i=0; i<left_by_gene.size()-1; i++)
			{
			myassert(pick <= all_infd);

			// split leftover into two groups
			all_infd -= left_by_gene[i];
			
			myassert(all_infd >= 0);

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
			if (!l->to->blocked)
				// add the stuff to the link and its end node
				l->to->frequencies[i] += add;
			}

		if (!l->to->blocked)
			// last allele gets the leftovers
			l->to->frequencies.back() += pick;
		// keep track
		left_all -= pick;
		left_by_gene.back() -= pick;
		}

	node->done = true;
	}


/** Run mechanistic genetics simulation on a range of nodes. */
template<class ITER, class BINOM_FUNC>
void annotate_frequencies_ibmm(const ITER & beg, const ITER & end, BINOM_FUNC & binom)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_frequencies_ibmm(*i, binom);

	for (ITER i=beg; i!=end; i++)
		(*i)->done = false;
	}


#endif	// IBMMIXED_H
