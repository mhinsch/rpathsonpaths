#ifndef NET_UTIL_H
#define NET_UTIL_H

#include <vector>

/** @file
 * Generic network utility code (doesn't require Rcpp or R). */

/** Find cycles in a network. */
struct Cycles
	{
	const vector<vector<size_t>> & net; //!< Network as children per node.
	vector<bool> visited;				//!< Within calls (cycles).
	vector<bool> done;					//!< Between calls (optimization).
	vector<size_t> stack;
	vector<vector<size_t>> res;			//! A list of cycles.

	/** Plain constructor. */
	Cycles(const vector<vector<size_t>> & network)
		: net(network), visited(network.size(), false), done(network.size(), false)
		{}

	/** Detect if there's at least one cycle in the subnetwork reachable from node cur. This 
	 * can be significantly faster than tracking down all cycles. Note that calling this
	 * directly on a non-source node might produce false positives down the line. */
	bool has_cycles(size_t cur)
		{
		visited[cur] = true;	// this will detect cycles
		done[cur] = true;		// we keep track of processed nodes so we can skip them

		// check children
		for (size_t i : net[cur])
			{
			if (visited[i])	// been here => cycle!
				return true;

			if (!done[i])	// hasn't been processed => do it now
				if (has_cycles(i))
					return true;
			}

		visited[cur] = false;	// entire subtree has been checked, clean up
		return false;
		}

	/** Find and record all cycles in the subnetwork reachable from cur. */
	void find_cycles(size_t cur)
		{
		stack.push_back(cur);	// keep track of current path
		done[cur] = true;

		// check children
		for (size_t i : net[cur])
			{
			// are we crossing our own path?
			const auto f = find(stack.begin(), stack.end(), i);
			// if so, add the entire loop to the list of cycles
			if (f != stack.end())
				{
				res.push_back(vector<size_t>(f, stack.end()));
				continue;
				}

			if (!done[i])
				find_cycles(i);
			}

		stack.pop_back();
		}
	};


#endif	// NET_UTIL_H
