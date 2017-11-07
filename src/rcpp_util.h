#ifndef RCPP_UTIL_H
#define RCPP_UTIL_H

#include <vector>
#include <unordered_map>
#include <string>
#include <set>

#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

#define R_ASSERT(cond, msg) {if (!(cond)) stop(msg);}

template<class T>
XPtr<T> make_S3XPtr(T * obj, const char * class_name, bool GC = true)
	{
	XPtr<T> xptr(obj, GC);
	xptr.attr("class") = class_name;
	return xptr;
	}


struct RRng
	{
	double outOf(double mi, double ma) const
		{
		return R::runif(mi, ma);
		}

	size_t operator()(size_t n) const
		{
		size_t r;

		while ((r = R::runif(0, n)) >= n);

		return r;
		}
	};


// collect a map of factor names to indices from an R factor
vector<size_t> adapt_factor(const IntegerVector & factor, vector<string> & names, 
	unordered_map<string, size_t> & idxs);


// wrap an R edge list in order to allow for somewhat sane handling of factors
class EdgeList
	{
	const IntegerVector & _from_raw, & _to_raw;
	vector<size_t> _from, _to;
	vector<string> _names;
	unordered_map<string, size_t> _idxs;

	size_t _f;

public:
	EdgeList(const IntegerVector & from, const IntegerVector & to)
		: _from_raw(from), _to_raw(to)
		{
		_f = int(from.inherits("factor")) + to.inherits("factor");

		R_ASSERT(_f != 1, "Both node lists have to be of the same type!");

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
