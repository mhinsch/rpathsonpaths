#ifndef IBMMIXED_H
#define IBMMIXED_H

#include <vector>

#include "util.h"

template<class NODE, class RNG>
void annotate_frequencies_ibmm(NODE * node, RNG & rng)
	{
	if (node->done)
		return;

	// we are pushing, so ignore leaves
	if (node->is_leaf() || node->rate_in <= 0)
		{
		node->normalize();
		node->done = true;
		return;
		}

	//std::cout << "ibm: " << node << "\n";

	for (auto link : node->inputs)
		annotate_frequencies_ibmm(link->from, rng);

	// no input set on this branch
	if (node->frequencies.size() == 0)
		{
		node->done = true;
		return;
		}

	// just in case, frequencies might have been initialized to a different scale than
	// input rate; we just assume input rate has priority
	if (node->is_root())
		node->normalize(node->rate_in_infd);

	const double outp = 
		std::accumulate(node->outputs.begin(), node->outputs.end(), 0.0, 
			[](const double & v, const auto & l){return v + l->rate;});

	if (outp <= 0)
		{
		node->normalize();
		node->done = true;
		return;
		}

	// need to use this since gross rates in the nodes will be inaccurate
	const double infd = 
		std::accumulate(node->frequencies.begin(), node->frequencies.end(), 0.0);

	if (infd <= 0)
		{
		node->done = true;
		return;
		}

	// reconstruct transmission rate
	// TODO this is ugly! *_ibmm should not depend on basic annotate_rates
	const double r_inf = 
		node->d_rate_in_infd / (node->rate_in - node->rate_in_infd + node->d_rate_in_infd);
	const double inp = node->rate_in;
	const double uninfd = inp - infd;
	// coin flip for each uninfected whether it becomes infected
	const int newly_infd = uninfd > 0 ? rng.binom(r_inf, int(uninfd)) : 0;

	myassert(r_inf >= 0 && inp > 0 && infd > 0 && uninfd >= 0 && newly_infd >=0);
	
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
			// if this is the last positice frequency
			const int add = rng.binom(std::min(1.0, p/rem), n);

			myassert(add >= 0);

			node->frequencies[i] += add;

			n -= add;
			rem -= p;

			myassert(n>=0 && rem>-0.0001); 
			}

		// R binom does weird stuff when rem and p are very close so we
		// skip the last step and just assign n directly
		node->frequencies.back() += n;
		}

// *** generate output
//
// We don't replace units that have been selected for output => we have to 
// use a multi-variate hypergeometric distribution (we actually use a regular 
// hypergeometric and the conditional method).

	for (auto l : node->outputs)
		l->to->frequencies.resize(node->frequencies.size(), 0.0);

	// we have to keep track of how many infected there are still left
	// all_infd == sum(node->frequencies)
	double all_infd = infd + newly_infd;
	double all_non_infd = inp - all_infd;
	auto left_by_gene = node->frequencies;

	//std::cout << "****\n";
	for (auto l : node->outputs)
		{
		l->rate_infd = 0;
		// how many go into this link (infd + non-infd)
		double pick = l->rate;
		// overall number of units left to pick from (for all links)
		double all_left = all_infd + all_non_infd;

		// std::cout << ">>" << all_left << " " << pick << "\n";
		// do alleles first, then after that assign rest to uninfected
		
		for (size_t i=0; i<left_by_gene.size(); i++)
			{
			// split leftover into two groups
			all_left -= left_by_gene[i];
			
			myassert(all_left >= 0);
			myassert(pick <= all_left+left_by_gene[i]);

			// picked units are either allele i or allele >i + uninfected
			const int add = rng.hypergeom(int(left_by_gene[i]), int(all_left), int(pick));
			
			//std::cout << left_by_gene[i] << " ";
			//std::cout << add << "\n";

			myassert(add >= 0);
			// adjust number to pick from for the next output
			left_by_gene[i] -= add;
			myassert(left_by_gene[i] >= 0);
			all_infd -= add;
			myassert(all_infd >= 0);
			// we pick less next time
			pick -= add;
			myassert(pick >= 0);
			// add the stuff to the link and its end node
			l->rate_infd += add;
			l->to->frequencies[i] += add;
			}
		// at this point all_left == all_non_infd
		// so, just assign the rest of the output's rate to non-infected
		all_non_infd -= pick;
		myassert(all_non_infd >= 0);
		}

	// adjust rates in the node; due to stochasticity these will be 
	// different from the gross rates calculated before
	node->rate_in_infd = infd + newly_infd;
	node->d_rate_in_infd = newly_infd;
	node->rate_out_infd = std::accumulate(node->outputs.begin(), node->outputs.end(), 0.0,
		[](auto s, const auto l){return s + l->rate_infd;});

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
