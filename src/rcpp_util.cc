#include "rcpp_util.h"


// collect a map of factor names to indices from an R factor
vector<size_t> adapt_factor(const IntegerVector & factor, vector<string> & names, 
	unordered_map<string, size_t> & idxs)
	{
	vector<size_t> nodes;

	const StringVector levels = factor.attr("levels");

	for (size_t f : factor)
		{
		if (IntegerVector::is_na(f))
			throw invalid_argument("missing value");

		// take into account 1-based indexing in R
		const string name = string(levels(f-1));

		// try inserting
		const auto ins_it = idxs.emplace(name, names.size());
		// success => name was new => new node
		if (ins_it.second)
			names.push_back(name);

		// Index of name can be != names.size() if the name already existed.
		nodes.push_back(ins_it.first->second);
		}

	// RVO
	return nodes;
	}


set<size_t> find_sinks(const EdgeList & el)
	{
	// all sources (i.e. all nodes in the from list
	vector<bool> is_source;

	// *** collect all sources

	for (size_t i=0; i<el.n_edges(); i++)
		{
		const size_t n = el.from(i);
		if (n >= is_source.size())
			// gaps are filled with false
			is_source.resize(n+1, false);
		// but this one is a source
		is_source[n] = true;
		}

	// *** now collect sinks

	set<size_t> scs;

	for (size_t i=0; i<el.n_edges(); i++)
		{
		const size_t n = el.to(i);
		// if it's not in the sources list is has to be a sink
		if (n >= is_source.size() || !is_source[n])
			scs.insert(n);
		}

	return scs;
	}

// works exactly like find_sinks but the other way around
set<size_t> find_sources(const EdgeList & el)
	{
	vector<bool> is_sink;

	for (size_t i=0; i<el.n_edges(); i++)
		{
		const size_t n = el.to(i);
		if (n >= is_sink.size())
			is_sink.resize(n+1, false);
		is_sink[n] = true;
		}

	set<size_t> scs;

	for (size_t i=0; i<el.n_edges(); i++)
		{
		const size_t n = el.from(i);
		if (n >= is_sink.size() || !is_sink[n])
			scs.insert(n);
		}

	return scs;
	}
