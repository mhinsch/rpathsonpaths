#include "rnet_util.h"

#include "libpathsonpaths/proportionalpick.h"

#include "rcpp_util.h"

void print_node_id(const Net_t * net, size_t i)
	{
	if (net->name_by_id.size())
		Rcout << net->name_by_id[i];
	else
		Rcout << i;
	}

void print_popsnode(const Node_t * n) {};


void _set_allele_freqs(Net_t * net, const List & ini)
	{
	const IntegerVector nodes = ini(0);
	const NumericMatrix freqs = ini(1);

	const size_t n_all = freqs.ncol();

	if (nodes.size() != freqs.nrow())
		stop("Invalid parameter 'iniDist': "
		"number of rows in $frequencies and number of elements in $nodes have to be equal!");
	
	for (auto n : net->nodes)
		n->frequencies.clear();

	const bool f = nodes.inherits("factor");

	StringVector levels = f ? nodes.attr("levels") : StringVector();

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = f ? net->id_by_name[string(levels(nodes(i)-1))] :
			nodes[i];
		if (n > net->nodes.size())
			stop("Invalid node index in iniDist$nodes!");

		net->nodes[n]->frequencies.resize(n_all, 0);

		for (size_t j=0; j<n_all; j++)
			net->nodes[n]->frequencies[j] = freqs(i, j);
		}
	}


size_t id_from_SEXP(const Net_t & net, SEXP id)
	{
	switch (TYPEOF(id))
		{
	case INTSXP:
		return as<int>(id);
	case STRSXP:
		return net.id_by_name.at(as<string>(id));
	default:
		stop("Node id has to be integer or string!");
		}
	}


void sample_node(const Node_t & node, size_t n, vector<size_t> & count)
	{
	if (count.size() != node.frequencies.size())
		stop("Invalid number of alleles in node!");

	ProportionalPick<> pick(0.000001, node.frequencies);
	RRng r;

	for (size_t i=0; i<n; i++)
		count[pick.pick(r)]++;
	}

