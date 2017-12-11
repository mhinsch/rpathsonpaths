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
 * different formats. If we are wrapping factors the main problem is that integer values
 * corresponding to names can differ between from and to. Therefore we have to re-index all
 * names and keep two internal lists of indices, plus a list of names and a map. For
 * non-factors we just keep the original integer vectors.
 *
 * @note EdgeList objects can *not* be kept across R function calls since the original R
 * vectors might get garbage collected in the meantime. */
class EdgeList
	{
	const IntegerVector & _from_raw, & _to_raw; //!< Raw data as obtained from R.
	/** Indices of the nodes' names (factor only). Needed since indices might differ in 
	 * the from and to factor.  */
	vector<size_t> _from, _to;					
	vector<string> _names;						//!<Union of the factors' levels (factor only).
	unordered_map<string, size_t> _idxs;		//!< Map from level to index (factor only).

	size_t _f;									//!< Whether we are wrapping a factor.

public:
	struct Edge
		{
		size_t from, to;
		Edge(size_t f, size_t t)
			: from(f), to(t)
			{}
		};

	/** Construct an EdgeList from two IntegerVector objects (which might be factors). */
	EdgeList(const IntegerVector & from, const IntegerVector & to)
		// we keep the raw data
		: _from_raw(from), _to_raw(to)
		{
		_f = int(from.inherits("factor")) + to.inherits("factor");

		R_ASSERT(_f != 1, "Both node lists have to be of the same type");
		R_ASSERT(from.size() == to.size(), "Not a valid edge list.");

		if (_f) // factors
			{
			// we use this to unify name indices in from and to factors (might differ).
			_from = adapt_factor(from, _names, _idxs);
			_to = adapt_factor(to, _names, _idxs);
			}
		else // not factors
			{
			// if the data is not factors we basically do nothing, but we still have to check
			// for NA
			for (size_t i=0; i<from.size(); i++)
				{
				R_ASSERT(!IntegerVector::is_na(from[i]), "missing value");
				R_ASSERT(!IntegerVector::is_na(to[i]), "missing value");
				}
			}
		}

	/** Get all names. Will be empty if not a factor. */
	vector<string> & names() 
		{
		return _names;
		}

	/** Get all names. Will be empty if not a factor. */
	const vector<string> & names() const
		{
		return _names;
		}

	/** Get map name -> index. Empty of not a factor */
	unordered_map<string, size_t> & idxs()
		{
		return _idxs;
		}

	/** Number of nodes. Slow if not a factor. */
	size_t n_nodes() const
		{
		if (_f)
			return _names.size();

		// find the highest index
		size_t n = 0;
		for (size_t i=0; i<_from_raw.size(); i++)
			{
			n = max(n, size_t(_from_raw[i]));
			n = max(n, size_t(_to_raw[i]));
			}

		return n+1;
		}

	/** Number of edges. */
	size_t n_edges() const
		{
		return _from_raw.size();
		}

	/** Is this a factor? */
	bool factor() const
		{
		return _f;
		}

	/** Index of ith from node. */
	size_t from(size_t i) const
		{
		return _f ? _from[i] : _from_raw(i);
		}

	/** Index of ith to node. */
	size_t to(size_t i) const
		{
		return _f ? _to[i] : _to_raw(i);
		}

	Edge edge(size_t i) const
		{
		return _f ? 
			Edge(_from[i], _to[i]) :
			Edge(_from_raw[i], _to_raw[i]);
		}

	/** Index of node by name. */
	size_t index(const string & name) const
		{
		return _idxs.at(name);
		}
	/** Name of node by index. */
	const string & name(size_t idx) const
		{
		return _names[idx];
		}

	/** Convert IntegerVector to factor if necessary. */
	void make_factor(IntegerVector & v) const
		{
		if (factor())
			{
			v.attr("class") = "factor";
			v.attr("levels") = names();
			}
		}
	
	struct EdgeIter 
		{
		typedef std::input_iterator_tag iterator_category;
		typedef Edge value_type;

		size_t i;
		const EdgeList & el;

		EdgeIter(size_t start, const EdgeList & edgelist)
			: i(start), el(edgelist)
			{}

		EdgeIter & operator++() {i++;}
		Edge operator*() const  {return el.edge(i);}
		bool operator!=(const EdgeIter & o) const {return o.i != i || &o.el != &el;}
		};	

	EdgeIter begin() const
		{
		return EdgeIter(0, *this);
		}
	EdgeIter end() const
		{
		return EdgeIter(_f ? _from.size() : _from_raw.size(), *this);
		}
	};

/** Find all sources in an edge list. */
set<size_t> find_sources(const EdgeList & el);
/** Find all sinks in an edge list. */
set<size_t> find_sinks(const EdgeList & el);

#endif	// RCPP_UTIL_H
