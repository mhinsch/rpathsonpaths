#ifndef RNETWORK_H
#define RNETWORK_H

/** @file Custom network class. */

#include "libpathsonpaths/transportnetwork.h"

#include <unordered_map>

using namespace std;

/** Our custom network class. Only necessary because we want to store factor stuff. */
template<class N, class L>
struct RNetwork : public TransportNetwork<N, L>
	{
	//! Map factor levels to internal node index.
	unordered_map<string, size_t> id_by_name;
	//! Factor level of each of our nodes.
	vector<string> name_by_id;
	};


#endif	// RNETWORK_H
