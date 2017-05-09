#include "rcpp_util.h"


// collect a map of factor names to indices from an R factor
vector<size_t> adapt_factor(const IntegerVector & factor, vector<string> & names, 
	unordered_map<string, size_t> & idxs)
	{
	vector<size_t> nodes;

	const StringVector levels = factor.attr("levels");

	for (size_t f : factor)
		{
		// take into account 1-based indexing in R
		const auto name = levels(f-1);

		// try inserting
		const auto ins_it = idxs.emplace(name, names.size());
		// success => name was new => new node
		if (ins_it.second)
			names.push_back(string(name));

		nodes.push_back(ins_it.first->second);
		}

	return nodes;
	}


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
