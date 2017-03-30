#ifndef PATHS_H
#define PATHS_H

#include <vector>

template<class NODE>
struct Path
	{
	typedef NODE node_t;

	vector<node_t *> seq;
	double p;

	Path(node_t * start = 0, double p_ini = 1.0)
		: seq(1, start), p(p_ini)
		{}

	size_t length() const
		{
		return seq.size();
		}

	const node_t * leaf() const
		{
		return seq.back();
		}
	};


template<class NODE>
struct PathNode : public NODE
	{
	typedef Path<PathNode<NODE> > path_t;
	vector<path_t> paths;
	};


template<class NODE>
void generate_paths(NODE * node)
	{
	typedef Path<NODE> path_t;
	// already done
	if (node->paths.size())
		return;

	vector<path_t> & paths = node->paths;

	// check all inputs
	for (const auto * link : node->inputs)
		{
		// sources start a path
		if (input->is_source())
			{
			// directly from source
			paths.push_back(path_t(node));
			// source -> colony
			paths.push_back(path_t(node, node->p_newly_infected()));
			// same node twice => colony
			paths.back().seq.push_back(node);
			}
		else 
			{
			NODE * from = link->from;
			
			// hasn't been done yet, process first
			if (from->paths.size() == 0)
				generate_paths(from);

			const double p_here = prob(from, node);

			// all paths ending at this input have a chance of going on to the 
			// current node
			for (const auto & p : from->paths)
				{
				// add this path to the current node
				paths.push_back(p);
				paths.back().seq.push_back(node);
				paths.back().p *= p_here;

				// and the colonization one as well
				paths.push_back(paths.back());
				// colonization is marked by repeated node
				paths.back().seq.push_back(node);
				paths.back().p *= node->p_newly_infected();
				}
			}

		}
	}

template<class ITER>
void generate_paths(const ITER & beg, const ITER & end)
	{
	for (ITER i=beg; i!=end; i++)
		generate_paths(*i);
	}
#endif	// PATHS_H

