#ifndef RNET_UTIL_H
#define RNET_UTIL_H

#include <vector>

#include <Rcpp.h>

#include "rpathsonpaths_types.h"

using namespace std;
using namespace Rcpp;

struct Drift
	{
	typedef typename Node_t::freq_t::value_type num_t;
	num_t theta;

	Drift(double t)
		: theta(t)
		{ }

	void operator()(const Node_t::freq_t & freqs, Node_t::freq_t & res)
		{
		if (res.size() != freqs.size())
			stop("Drift: result vector has to be same size as input vector!");

		num_t norm = 0.0;		

		for (size_t i=0; i<freqs.size(); i++)
			norm += (res[i] = R::rgamma(freqs[i] * theta, 1.0));

		for (num_t & f : res)
			f /= norm;
		}
	};

// binomial and hypergeometric dists for the ibm
struct Rng
	{
	double binom(double p, double n) const
		{
		return R::rbinom(n, p);
		}

	double hypergeom(double n1, double n2, double k) const
		{
		return R::rhyper(n1, n2, k);
		}
	};

void print_node_id(const Net_t * net, size_t i);
void print_popsnode(const Node_t * n);

void _set_allele_freqs(Net_t * net, const List & ini);

size_t id_from_SEXP(const Net_t & net, SEXP id);

void sample_node(const Node_t & node, size_t n, vector<size_t> & count);



#endif	// RNET_UTIL_H
