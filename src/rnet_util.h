#ifndef RNET_UTIL_H
#define RNET_UTIL_H

#include <vector>

#include <Rcpp.h>

#include "libpathsonpaths/proportionalpick.h"

#include "rpathsonpaths_types.h"
#include "rcpp_util.h"

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
	int binom(double p, int n) const
		{
		return R::rbinom(n, p);
		}

	int hypergeom(int n1, int n2, int k) const
		{
		return R::rhyper(n1, n2, k);
		}
	};

void print_node_id(const Net_t * net, size_t i);
void print_popsnode(const Node_t * n);

void _set_allele_freqs(Net_t * net, const List & ini);

size_t id_from_SEXP(const Net_t & net, SEXP id);

void sample_node(const Node_t & node, size_t n, vector<size_t> & count);

template<class CONT>
void sample_alleles_node(const Node_t & node, CONT & alleles)
	{
	ProportionalPick<> pick(0.000001, node.frequencies);
	RRng r;

	for (auto & a : alleles)
		a = pick.pick(r);
	}


double distance_SNP(const Node_t & n1, const Node_t & n2);
double distance_freq(const Node_t & n1, const Node_t & n2);
double distance_EHamming(const Node_t & n1, const Node_t & n2);

#endif	// RNET_UTIL_H
