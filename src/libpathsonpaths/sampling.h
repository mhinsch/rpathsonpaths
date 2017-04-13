#ifndef SAMPLING_H
#define SAMPLING_H

#include "transportgraph.h"
#include "proportionalpick.h"

template<class RNG>
void pick_sample(const Node * node, RNG & rng, vector<Path *> & res, size_t n_samples)
	{
	ProportionalPick pp(0.0001);

	pp.setup(begin(node.paths), end(node.paths), 
		[](decltype(begin(node.paths)) p_iter){return p_iter->p;});

	for (size_t i=0; i<n_samples; i++)
		{
		res.push_back(&node->paths[pp.pick(rng)]);
		}
	}

struct PhyloNode
	{
	vector<Path*> paths;
	Node * node;
	PhyloNode * ancestor;
	vector<PhyloNode*> children;
	size_t dist_from_root;

	PhyloNode()
		: paths(), node(0), ancestor(0), children(0), dist_from_root(0)
		{}
	};


struct Joint
	{
	const Path * joiner;
	vector<const Path *> joinees;
	size_t time;

	Joint(Path * j, Path * je, size_t t)
		: joiner(j), time(t)
		{
		joinees.push_back(je);
		}
	};


void build_phylogeny(const vector<Path *> & paths)
	{
	const size_t n_paths = size(paths);

	vector<Joint *> joints;
	
	vector<vector<Joint *> > joining;

	joining.resize(n_paths, 0);

	for (size_t i=0; i<npaths-1; i++)
		{
		const Path * pi = paths[i];

		for (size_t j=i+1; j<npaths; j++)
			{
			const Path * pj = paths[j];

			const size_t ovl = count_common_root(pi, pj);
			// completely distinct
			if (ovl == 0)
				continue;

			// pi has a colonization event at this point
			if (pi->seq[ovl] == pi->seq[ovl-1])
				{
				// check all of i's colon. events to see if we have this one already
				auto jt = find_if(joining[i].begin(), joining[i].end(), 
					[ovl](Joint * joint){return joint->time == ovl});

				if (jt == joining[i].end())
					{
					joints.push_back(new Joint(pi, pj, ovl));
					joining[i].push_pack(joints.back());
					}
				else
					(*jt)->joinees.push_back(pj);
				}

			// pj has a col. event 
			if (pj->seq[ovl] == pj->seq[ovl-1])
				{
				// check all of j's colon. events
				auto jt = find_if(joining[j].begin(), joining[j].end(), 
					[ovl](Joint * joint){return joint->time == ovl});

				if (jt == joining[j].end())
					{
					joints.push_back(new Joint(pj, pi, ovl));
					joining[j].push_pack(joints.back());
					}
				else
					(*jt)->joinees.push_back(pi);
				}
			}
		}

	sort(joints.begin(), joints.end(), [](Joint * j1, Joint * j2){return j1->time > j2->time});

	for (Joint * joint : joints)
		{
		// TODO where does p_join come from?
		const double p_none = pow(1-p_join, joint->joinees.size());
		
		// TODO where does rng come from?
		// actually joining
		if (!rng.choice(p_none))
			{
			size_t joinee = rng(joint->joinees.size());

			
			}
		}
	}

// find common start sequence between two paths
// returns position of first unequal element
// TODO make generic for all types of sequences
template<class SEQ>
size_t count_common_root(const SEQ & s1, const SEQ & sp2)
	{
	const size_t st = min(size(s1), size(s2));

	auto i1 = begin(s1);
	auto i2 = begin(s2);

	for (size_t i=0; i<st; i++, i1++, i2++)
		if (*i1 != *i2)
			return i;

	return st;
	}


template<class ELEM, class SEQ>
void remove(ELEM & el, SEQ & seq)
	{
	swap(el, seq.back());
	seq.pop_back();
	}

#endif	// SAMPLING_H
