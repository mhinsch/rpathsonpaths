#ifndef NET_UTIL_H
#define NET_UTIL_H

#include <vector>

struct Cycles
	{
	const vector<vector<size_t>> & net;
	vector<bool> visited;
	vector<bool> done;
	vector<size_t> stack;
	vector<vector<size_t>> res;

	Cycles(const vector<vector<size_t>> & network)
		: net(network), visited(network.size(), false), done(network.size(), false)
		{}

	bool has_cycles(size_t cur)
		{
		visited[cur] = true;
		done[cur] = true;

		for (size_t i : net[cur])
			{
			if (visited[i])
				return true;

			if (!done[i])
				if (has_cycles(i))
					return true;
			}

		visited[cur] = false;
		return false;
		}

	void find_cycles(size_t cur)
		{
		stack.push_back(cur);
		done[cur] = true;

		for (size_t i : net[cur])
			{
			const auto f = find(stack.begin(), stack.end(), i);
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
