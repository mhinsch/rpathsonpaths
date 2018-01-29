#ifndef NET_UTIL_H
#define NET_UTIL_H

#include <vector>
#include <algorithm>
#include <unordered_set>

using std::vector;

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
			const auto f = std::find(stack.begin(), stack.end(), i);
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


/** Generate a random scale-free network. The function uses the Barabasi-Albert 
 preferential attachment algorithm, slightly modified to allow for directedness and 
 isolated initial nodes.

 @param from, to Containers that the generated edge list will be written to.
 @param n_nodes Number of (non-source) nodes to generate.
 @param n_sources Number of source nodes to initialize the network with (has to
 be at least 1). Note that there is no guarantee all source nodes will become
 part of the network.  
 @param m_dist The probability distribution to draw the number of inputs for
 new nodes from. m_dist has to be a function object that receives the node
 index and returns the number of nodes.
 @param zero_appeal Constant to be added to the nodes' attractiveness.
 @param rng A random number generator.
 @param compact Whether to remove isolated source nodes. */
template<class INT_CONT, class DIST, class RNG>
void net_gen_prefattach(INT_CONT & from, INT_CONT & to, int n_nodes, int n_sources,
	const DIST & m_dist, float zero_appeal, RNG & rng, bool compact=false)
	{
	vector<float> weight(n_nodes + n_sources, 0);

	from.reserve(n_nodes);
	to.reserve(n_nodes);
	
	for (size_t i=0; i<n_sources; i++)
		weight[i] = zero_appeal;

	// we need to keep track of the sum of all weights
	double sum = n_sources * zero_appeal;


	for (size_t i=0; i<n_nodes; i++)
		{
		// index of current node
		const size_t node = i + n_sources;
		// how many inputs
		const size_t n_inp = m_dist(i) + 1;

		for (size_t j=0; j<n_inp; j++)
			{
			// random node
			size_t r_inp = rng.outOf(0, sum);
			size_t inp = 0;
			// find previous node in weight list
			while (r_inp > weight[inp])
				r_inp -= weight[inp++];

			from.push_back(inp);
			to.push_back(node);

			// input node gains a connection
			weight[inp]++;
			// and sum increases accordingly
			sum++;
			}

		// new node has 0 outputs
		weight[node] = zero_appeal;
		sum++;
		}


	if (compact)
		{
		// *** remove isolated nodes === make node indices contiguous
		//     we re-use weight to store how much we have to count the index for a given
		//     node down by

		int reduce = 0;

		// find isolated sources
		for (size_t i=0; i<n_sources; i++)
			{
			// sources that are still at zero_appeal have no outputs => isolated
			if (weight[i] == zero_appeal)
				reduce++;
			else
				weight[i] = reduce;
			}

		// regular nodes can't be isolated, but we still have to change their index 
		fill(weight.begin()+n_sources, weight.end(), reduce);

		for (size_t i=0; i<from.size(); i++)
			{
			from[i] -= weight[from[i]];
			to[i] -= weight[to[i]];
			}
		}
	}


/** Identify separate sub-networks in a network described as an edge list.

   @param beg, end Iterators that point to the beginning and end of an edge list.

   @return An integer vector with the sub-network id of each node. Note that id's start at
   1 and are not guaranteed to be contiguous. */
template<class EI>
vector<int> colour_network(EI beg, EI end)
	{
	// colour of nodes
	vector<int> colour;

	int next_col = 1;

	for (; beg != end; ++beg)
		{
		const size_t f = (*beg).from, t = (*beg).to;

		if (max(f, t) >= colour.size())
			colour.resize(max(f, t)+1, 0);

		if (colour[f] == colour[t])
			{
			// not coloured yet, colour them
			if (colour[f] == 0)
				{
				colour[f] = next_col++;
				colour[t] = colour[f];
				}
			// otherwise they are the same colour which is also fine
			}
		else
			{
			// one of them is not coloured => take the other one's colour
			if (colour[f] == 0)
				{
				colour[f] = colour[t];
				continue;
				}
			if (colour[t] == 0)
				{
				colour[t] = colour[f];
				continue;
				}
			// two different colours, have to change all instances of one of them
			const int oldc = colour[t];
			const int newc = colour[f];

			replace(colour.begin(), colour.end(), oldc, newc);
			}
		}

	return colour;
	}


template<class D>
extern void set_value(D & d, size_t x, size_t y, int v);

/** Determine pairwise topological distances on a list of nodes.
 * @tparam CONT Node container type.
 * @tparam DIST Dist matrix type.
 * @param nodes Container of pointers to nodes. If nodes is sorted the function will 
 * run slightly faster.
 * @param dists Matrix of distances. Calls element & at(DIST &, size_t, size_t) for 
 * element access.
 * @param leaves_only Save time by generating only distances between leaf nodes. The rest
 * of the distance matrix is filled with -1 or 0 (diagonal).
 */
template<class CONT, class DIST>
void distances(const CONT & nodes, DIST & dists, bool leaves_only = true)
	{
	typedef typename CONT::value_type NODEP;

	const bool sorted = is_sorted(nodes.begin(), nodes.end());

	for (size_t i=0; i<nodes.size(); i++)
		for (size_t j=0; j<nodes.size(); j++)
			set_value(dists, i, j, i==j ? 0 : -1);

	unordered_set<NODEP> visited;
	vector<NODEP> stack_next, stack_cur;

	// find all distances for each of the leaf nodes
	for (const auto start : nodes)
		{
		size_t dist = 1;

		// new start node, thus clear everything
		visited.clear();
		stack_cur.clear();
		// first one to process		
		stack_cur.push_back(start);

		const size_t idx_start = (sorted ? 
				std::lower_bound(nodes.begin(), nodes.end(), start) :
				std::find(nodes.begin(), nodes.end(), start))
			- nodes.begin();

		while (stack_cur.size())
			{
			for (const auto n : stack_cur)
				{
				// there might be cycles, so this is possible
				if (visited.count(n))
					continue;

				for (const auto link : n->inputs)
					{
					NODEP parent = link->from;
					// parent has been visited, skip
					if (visited.count(parent))
						continue;

					if (!leaves_only)
						{
						const size_t idx_parent = (sorted ?
								std::lower_bound(nodes.begin(), nodes.end(), parent) :
								std::find(nodes.begin(), nodes.end(), parent))
							- nodes.begin();

						set_value(dists, idx_start, idx_parent, dist);
						set_value(dists, idx_parent, idx_start, dist);
						}

					stack_next.push_back(parent);
					}

				for (const auto link : n->outputs)
					{
					NODEP child = link->to;
					// child has been visited, skip
					if (visited.count(child))
						continue;

					if (!leaves_only || child->is_leaf())
						{
						const size_t idx_child = (sorted ?
								std::lower_bound(nodes.begin(), nodes.end(), child) :
								std::find(nodes.begin(), nodes.end(), child))
							- nodes.begin();

						set_value(dists, idx_start, idx_child, dist);
						set_value(dists, idx_child, idx_start, dist);
						}

					stack_next.push_back(child);
					}

				// we assume there are no links to self, so it should
				// be fine doing this here
				visited.insert(n);
				}

			// next stack becomes current
			stack_cur.swap(stack_next);
			// and clear the next one so that we can fill it again
			stack_next.clear();
			// one layer out
			dist++;
			}
		}
	}
#endif	// NET_UTIL_H
