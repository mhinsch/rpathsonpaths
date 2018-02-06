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

		if (std::max(f, t) >= colour.size())
			colour.resize(std::max(f, t)+1, 0);

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


template<class INP>
struct ret_type
	{
	typedef double T;
	};

// slightly convoluted, but the only way I've found to make this work for 
// NumericMatrix *and* boost::multi_array
template<class D>
typename ret_type<D>::T & at(D & d, size_t x, size_t y);

/** Determine pairwise topological distances on a list of nodes.
 * @tparam CONT Node container type.
 * @tparam DIST Dist matrix type.
 * @param nodes Container of pointers to nodes. If nodes is sorted the function will 
 * run slightly faster.
 * @param dists Matrix of distances. Calls element & at(DIST &, size_t, size_t) for 
 * element access.
 */
template<class CONT, class DIST>
void distances(const CONT & nodes, DIST & dists)
	{
	// In principle we go through all nodes and for each of them do a full width-first 
	// search of the network. However, we know that for two nodes A and B and a third node
	// C the shortest path from A to B either passes through C or it doesn't. If it does
	// and we know all shortest paths from C as well as the shortest path from A to C, we know
	// the length of A to B. If it doesn't we will encounter it at some other point in the 
	// search. In either case, however, there is no need to search beyond C.

	typedef typename CONT::value_type NODEP;

	const bool sorted = is_sorted(nodes.begin(), nodes.end());

	auto get_node_idx = [&nodes, sorted](const NODEP n)
		{
		if (sorted)
			{
			const auto it = std::lower_bound(nodes.begin(), nodes.end(), n);
			return (it == nodes.end() || *it != n ? nodes.end() : it) - nodes.begin();
			}

		return std::find(nodes.begin(), nodes.end(), n) - nodes.begin();
		};

	static std::unordered_set<NODEP> skip;
	static std::unordered_set<NODEP> done_NIM;
	static vector<NODEP> stack_next, stack_cur;
	stack_next.clear();

	auto process_node = [&] (const NODEP node, size_t idx_start, int dist)
		{
		// only needed here because it won't ever be stacked
		if (skip.count(node))
			return;

		const size_t idx_node = get_node_idx(node);
		// has been processed before
		// NOTE: never triggers for nodes not included in nodes list 
		// (e.g. non-leaves if only leaves are processed)
		if (idx_start < nodes.size() && idx_node < idx_start)
			{
			for (size_t j=0; j<nodes.size(); j++)
				{
				const auto d = at(dists, idx_node, j);
				if (d > 0)
					{
					auto & e = at(dists, idx_start, j);
					if (e < 0 || e > d+dist)
						{
						e = d+dist;
						at(dists, j, idx_start) = d+dist;
						}
					}
				}

			// ignore from now on
			skip.insert(node);
			return;
			}

		// node not in list => we can't use matrix to check if we did it already
		if (idx_node == nodes.size())
			{
			if (!done_NIM.count(node))
				{
				stack_next.push_back(node);
				done_NIM.insert(node);
				}
			return;
			}

		auto & v = at(dists, idx_start, idx_node);

		// either not set or higher => set and process next turn
		if (v < 0 || v > dist)
			{
			v = dist;
			at(dists, idx_node, idx_start) = dist;
			stack_next.push_back(node);
			}
		};

	for (size_t i=0; i<nodes.size(); i++)
		for (size_t j=0; j<nodes.size(); j++)
			at(dists, i, j) = i==j ? 0 : -1;

	// find all distances for each of the nodes
	for (auto i=nodes.begin(); i!=nodes.end(); i++)
		{
		size_t dist = 1;

		skip.clear();
		done_NIM.clear();
		// new start node, thus clear everything
		stack_cur.clear();
		// first one to process		
		stack_cur.push_back(*i);
		const size_t idx_start = i - nodes.begin();
		// not in matrix => remember this one is done
		if (i == nodes.end())
			done_NIM.insert(*i);

		while (stack_cur.size())
			{
			for (const auto n : stack_cur)
				{
				for (const auto link : n->inputs)
					{
					const NODEP parent = link->from;
					process_node(parent, idx_start, dist);
					}

				for (const auto link : n->outputs)
					{
					const NODEP child = link->to;
					process_node(child, idx_start, dist);
					}
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
