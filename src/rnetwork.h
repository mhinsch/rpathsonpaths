#ifndef RNETWORK_H
#define RNETWORK_H

#include "libpathsonpaths/network.h"

#include <unordered_map>

using namespace std;

template<class N, class L>
struct RNetwork : public Network<N, L>
	{
	//! map factor levels to internal node index
	unordered_map<string, size_t> id_by_name;
	//! factor level of each of our nodes
	vector<string> name_by_id;
	};


#endif	// RNETWORK_H
