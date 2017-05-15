#ifndef TRANSPORTNETWORK_H
#define TRANSPORTNETWORK_H

#include "network.h"

template<class N, class L>
struct TransportNetwork : public Network<N, L>
	{
	/** Make node @a s an external source with proportion of infected set to @a prop_infd. */
	void set_source(size_t s, double prop_infd)
		{
		assert(this->nodes.size() > s && this->nodes[s] != 0);
		
		this->nodes[s]->rate_in = 1.0;
		this->nodes[s]->rate_in_infd = prop_infd;
		this->nodes[s]->d_rate_in_infd = 0;
		}
	};

#endif	// TRANSPORTNETWORK_H
