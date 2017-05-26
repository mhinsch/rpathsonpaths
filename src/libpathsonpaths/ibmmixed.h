#ifndef IBMMIXED_H
#define IBMMIXED_H

#include <vector>


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

	for (auto * link : node->inputs)
		annotate_frequencies_ibmm(link->from, rng);

	// no input set on this branch
	if (node->frequencies.size() == 0)
		{
		node->done = true;
		return;
		}

	const double outp = 
		std::accumulate(node->outputs.begin(), node->outputs.end(), 0.0, 
			[](const double & v, const auto & l){return v + l->rate;});

	if (outp <= 0)
		{
		node->normalize();
		node->done = true;
		return;
		}

	// reconstruct transmission rate
	// TODO this is ugly! *_ibmm should not depend on basic annotate_rates
	const double r_inf = 
		node->d_rate_in_infd / (node->rate_in - node->rate_in_infd + node->d_rate_in_infd);
	const double inp = node->rate_in;
	// need to use this since gross rates in the nodes will be inaccurate
	const double infd = 
		std::accumulate(node->frequencies.begin(), node->frequencies.end(), 0.0);
	const double uninfd = inp - infd;
	// coin flip for each uninfected whether it becomes infected
	const double newly_infd = uninfd > 0 ? rng.binom(r_inf, uninfd) : 0;
	
// *** transmission 
//
// Units that transmitted infection remain in the pool => polynomial distribution
// (or binomial + conditional method).

	// if there are newly infected units, distribute them according to prop. of 
	// alleles
	// multinomial assumes that all infections happen simultaneously
	if (newly_infd > 0)
		{
		double n = newly_infd;
		double rem = 1.0;
		for (size_t i=0; i<node->frequencies.size() && n>0 && rem>0; i++)
			{
			const double p = node->frequencies[i] / infd;
			const double add = rng.binom(p/rem, n);

			node->frequencies[i] += add;

			n -= add;
			rem -= p;
			}
		}

// *** generate output
//
// We don't replace units that have been selected for output => we have to 
// use a multi-variate hypergeometric distribution (we actually use a regular 
// hypergeometric and the conditional method).

	for (auto * l : node->outputs)
		l->to->frequencies.resize(node->frequencies.size());

	// we have to keep track of how many infected there are still left
	// all_infd == sum(node->frequencies)
	double all_infd = infd + newly_infd;
	double all_non_infd = inp - all_infd;
	auto left_by_gene = node->frequencies;

	for (auto * l : node->outputs)
		{
		l->rate_infd = 0;
		// how many go into this link
		double pick = l->rate;
		// overall number of units left to pick from
		double all_left = all_infd + all_non_infd;

		// do alleles first, then after that assign rest to uninfected
		
		for (size_t i=0; i<left_by_gene.size(); i++)
			{
			// split leftover into two groups
			all_left -= left_by_gene[i];
			// picked units are either allele i or allele >i + uninfected
			int add = rng.hypergeom(left_by_gene[i], all_left, pick);
			// adjust number to pick from for the next output
			left_by_gene[i] -= add;
			all_infd -= add;
			// we pick less next time
			pick -= add;
			// add the stuff to the link and its end node
			l->rate_infd += add;
			l->to->frequencies[i] += add;
			}
		// at this point all_left == all_non_infd
		// so, just assign the rest of the output's rate to non-infected
		all_non_infd -= pick;
		}

	// adjust rates in the node; due to stochasticity these will be 
	// different from the gross rates calculated before
	node->rate_in_infd = infd + newly_infd;
	node->d_rate_in_infd = newly_infd;
	node->rate_out_infd = 0;
	for (const auto * l : node->outputs)
		node->rate_out_infd += l->rate_infd;

	// TODO normalize frequencies
	node->done = true;
	}


template<class ITER, class BINOM_FUNC>
void annotate_frequencies_ibmm(const ITER & beg, const ITER & end, BINOM_FUNC & binom)
	{
	for (ITER i=beg; i!=end; i++)
		annotate_frequencies_ibmm(*i, binom);
	}


#endif	// IBMMIXED_H
