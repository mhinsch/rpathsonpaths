#ifndef RNET_UTIL_H
#define RNET_UTIL_H

/** @file Network utilities that require R/Rcpp. */


#include <vector>

#include <Rcpp.h>

#include "libpathsonpaths/proportionalpick.h"

#include "rpathsonpaths_types.h"
#include "rcpp_util.h"
#include "net_util.h"

using namespace std;
using namespace Rcpp;


/** Drift operator for the continuous model. Implements a Dirichlet distribution using
 * R::rgamma. */
struct Drift
	{
	typedef typename Node_t::freq_t::value_type num_t;
	num_t theta;	//!< Scaling parameter for the Dirichlet distribution.

	Drift(double t)
		: theta(t)
		{ }

	/** Apply drift to dreqs and store result in res. */
	void operator()(const Node_t::freq_t & freqs, Node_t::freq_t & res)
		{
		R_ASSERT(res.size() == freqs.size(), 
			"Drift: result vector has to be same size as input vector");

		num_t norm = 0.0;		

		// draw from a Gamma distribution
		for (size_t i=0; i<freqs.size(); i++)
			norm += (res[i] = R::rgamma(freqs[i] * theta, 1.0));

		// normalize
		for (num_t & f : res)
			f /= norm;
		}
	};


/** Binomial and hypergeometric distributions for the mechanistic model. */
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


/** Print (using R output) node i of network net. */
void print_node_id(const Net_t * net, size_t i);

/** Print a node. Does nothing at the moment. */
void print_popsnode(const Node_t * n);


/** Set a bunch of allele frequencies. ini has to contain a list of node ids as first and a
 * matrix (node x alleles) of allele frequencies as second element. */
void _set_allele_freqs(Net_t * net, const List & ini);


/** Get a node id from an R SEXP containing either an integer or a string (for factors). */
size_t id_from_SEXP(const Net_t & net, SEXP id);


/** Obtain a number of random samples from a node, storing the number of times each allele
 * was drawn in count. */
void sample_node(const Node_t & node, size_t n, vector<size_t> & count);

/** Obtain a number of random samples from a node, storing the allele id of each draw in 
 * alleles. */
template<class CONT>
void sample_alleles_node(const Node_t & node, CONT & alleles)
	{
	ProportionalPick<> pick(0.000001, node.frequencies);
	RRng r;

	for (auto & a : alleles)
		a = pick.pick(r);
	}


/** Mean square difference in allele frequencies between two nodes. */
double distance_freq(const Node_t & n1, const Node_t & n2);

/** Expected value of the Hamming distance between two nodes. */
double distance_EHamming(const Node_t & n1, const Node_t & n2);

/** Template specialization for NumericMatrix. */
template<>
//inline auto at(NumericMatrix & m, size_t x, size_t y) -> decltype(*m.begin()) &
inline double & at(NumericMatrix & m, size_t x, size_t y)
	{
	return m.at(x, y);
	}

#endif	// RNET_UTIL_H
