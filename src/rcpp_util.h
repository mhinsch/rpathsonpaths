#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

/** @file
 * Some util code that requires Rcpp.
 */

#include <vector>
#include <unordered_map>
#include <string>
#include <set>

#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

/** Assert that forwards to R's stop. i
 * @param cond Condition to evaluate. 
 * @param msg Error message. */
#define R_ASSERT(cond, msg) ((cond) ? (void)0 : stop(msg))

/** Wrap an object into an R S3 pointer.
 * @tparam Class of the object.
 * @param obj Object to wrap.
 * @param class_name S3 class name.
 * @param GC Whether to garbage collect the object (Rcpp internal). */
template<class T>
XPtr<T> make_S3XPtr(T * obj, const char * class_name, bool GC = true)
	{
	XPtr<T> xptr(obj, GC);
	xptr.attr("class") = class_name;
	return xptr;
	}

/** Shallow wrapper for the R RNG. */
struct RRng
	{
	/** Generate a uniform random double. */
	double outOf(double mi, double ma) const
		{
		return R::runif(mi, ma);
		}

	/** Generate a uniform random integer. */
	size_t operator()(size_t n) const
		{
		size_t r;

		while ((r = R::runif(0, n)) >= n);

		return r;
		}
	};


/** Collect a map of factor names to indices from an R factor.
 * @param factor An R factor to adapt.
 * @param names Will contain all names in the factor.
 * @param idxs Will contain a map from names to indices.
 * @return An ordered list of indices. */
vector<size_t> adapt_factor(const IntegerVector & factor, vector<string> & names, 
	unordered_map<string, size_t> & idxs);


/** A pragmatic wrapper for edgelists in either factor or integer vector format. 
 *
 * This wrapper class hides most of the complexities of having to deal with edgelists in two
 * different formats. */
class EdgeList
	{
	const IntegerVector & _from_raw, & _to_raw; //!< Raw data as obtained from R.
	vector<size_t> _from, _to;					//!< 
	vector<string> _names;						//!< The factor's levels.
	unordered_map<string, size_t> _idxs;		//!< Map from level to index.

	size_t _f;									//!< Whether we are wrapping a factor.

public:
	EdgeList(const IntegerVector & from, const IntegerVector & to)
		: _from_raw(from), _to_raw(to)
		{
		_f = int(from.inherits("factor")) + to.inherits("factor");

		R_ASSERT(_f != 1, "Both node lists have to be of the same type");

		if (_f)
			{
			_from = adapt_factor(from, _names, _idxs);
			_to = adapt_factor(to, _names, _idxs);
			}
		else
			{
			for (size_t i=0; i<from.size(); i++)
				{
				R_ASSERT(!IntegerVector::is_na(from[i]), "missing value");
				R_ASSERT(!IntegerVector::is_na(to[i]), "missing value");
				}
			}
		}

	vector<string> & names() 
		{
		return _names;
		}

	unordered_map<string, size_t> & idxs()
		{
		return _idxs;
		}

	size_t n_nodes() const
		{
		if (_f)
			return _names.size();

		// strictly speaking this gives us the highest index
		size_t n = 0;
		for (size_t i=0; i<_from_raw.size(); i++)
			{
			n = max(n, size_t(_from_raw[i]));
			n = max(n, size_t(_to_raw[i]));
			}

		return n+1;
		}

	size_t n_edges() const
		{
		return _from_raw.size();
		}

	bool factor() const
		{
		return _f;
		}

	size_t from(size_t i) const
		{
		return _f ? _from[i] : _from_raw(i);
		}

	size_t to(size_t i) const
		{
		return _f ? _to[i] : _to_raw(i);
		}

	size_t index(const string & name) const
		{
		return _idxs.at(name);
		}
	const string & name(size_t idx) const
		{
		return _names[idx];
		}
	};


set<size_t> find_sources(const EdgeList & el);
set<size_t> find_sinks(const EdgeList & el);


#endif	// RCPP_UTIL_H
