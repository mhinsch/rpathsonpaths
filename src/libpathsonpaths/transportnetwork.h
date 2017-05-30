#ifndef TRANSPORTNETWORK_H
#define TRANSPORTNETWORK_H

#include "util.h"

#include "network.h"


template<class N, class L>
struct TransportNetwork : public Network<N, L>
	{
	/** Make node @a s an external source with rate of infected set to @a r_infd. */
	void set_source(size_t s, double r_infd, double r_in = 1.0)
		{
		myassert(this->nodes.size() > s && this->nodes[s] != 0);
		
		this->nodes[s]->rate_in = r_in;
		this->nodes[s]->rate_in_infd = r_infd;
		this->nodes[s]->d_rate_in_infd = 0;
		}
	};

#endif	// TRANSPORTNETWORK_H
