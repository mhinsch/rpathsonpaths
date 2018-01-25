#include "rnet_util.h"

#include "libpathsonpaths/sputil.h"


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

	R_ASSERT(nodes.size() == freqs.nrow(), "Invalid parameter 'iniDist': "
		"number of rows in frequencies and number of elements in nodes have to be equal");	

	// init and reset nodes
	for (auto n : net->nodes)
		{
		// set all to 0
		n->frequencies.resize(n_all);
		fill(n->frequencies.begin(), n->frequencies.end(), 0);
		// root nodes start with wild type
		if (n->is_root())
			n->frequencies[0] = 1.0;
		}

	const bool f = nodes.inherits("factor");

	StringVector levels = f ? nodes.attr("levels") : StringVector();

	// assign ini frequencies
	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = f ? net->id_by_name[string(levels(nodes(i)-1))] : nodes[i];

		R_ASSERT(n < net->nodes.size(), "Invalid node index");

		Node_t * node = net->nodes[n];

		for (size_t j=0; j<n_all; j++)
			node->frequencies[j] = freqs(i, j);

		// no additional input into this node
		node->blocked = true;
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
		stop("Node id has to be integer or string");
		}
	}


void sample_node(const Node_t & node, size_t n, vector<size_t> & count)
	{
	R_ASSERT(count.size() == node.frequencies.size(), "Invalid number of alleles in node");

	ProportionalPick<> pick(0.000001, node.frequencies);
	RRng r;

	for (size_t i=0; i<n; i++)
		count[pick.pick(r)]++;
	}


double distance_freq(const Node_t & n1, const Node_t & n2)
	{
	double d = 0.0;

	for (int i=0; i<n1.frequencies.size(); i++)
		d += pow<2>(n1.frequencies[i] - n2.frequencies[i]);

	return d/n1.frequencies.size();
	}


double distance_EHamming(const Node_t & n1, const Node_t & n2)
	{
	double d = 0.0;

	for (int i=0; i<n1.frequencies.size(); i++)
		d += n1.frequencies[i] * n2.frequencies[i];

	return 1.0 - d;
	}


template<>
double & at(NumericMatrix & m, size_t x, size_t y)
	{
	return m(x, y);
	}

template<>
double at(const NumericMatrix & m, size_t x, size_t y)
	{
	return m(x, y);
	}
