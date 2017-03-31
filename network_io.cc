#include "network_io.h"

#include "sputil.h"

using namespace std;

void read_network(istream & in, AbstractNetwork & net)
	{
	istringstream istr;
	string str;

	do
		{
		skip_space(in, str);
		if (!in.good())
			return;
		
		istr.clear();
		istr.str(str);
		string tag;
		size_t from, to;
		double r;

		istr >> tag;

		istr >> from;
		istr >> to;
		istr >> r;

		if (tag == "N")
			net.add_link(from, to, r);
		else if (tag == "S")
			net.set_source(from, r);
		else
			{
			cerr << "Error: unknown node type " <<tag << "!\n";
			exit(1);
			}
		}
	while(in.good());
	}

// TODO fix
/*void write_network(ostream & out, const AbstractNetwork & net)
	{
	for (size_t i=0; i<net.n_nodes(); i++)
		{
		const Network::node_t  & inps = net.node(i);

		for (size_t j=0; j<inps.size(); j++)
			{
			const Input & inp = inps[j];

			out << (inp.source ? "S" : "N") << "\t";
			out << inp.node << "\t";
			out << i << "\t";
			out << inp.rate << "\n";
			}	
		}
	} */
