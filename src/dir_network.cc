#include "dir_network.h"

#include "libpathsonpaths/proportionalpick.h"

XPtr<Net_t> PopsNetwork(DataFrame links, DataFrame external, double transmission)
	{
	const IntegerVector inputs = links["inputs"];
	const NumericVector rates = links["rates"];
	const IntegerVector outputs = links["outputs"];

	const IntegerVector ext = external["nodes"];
	const NumericVector ext_rates = external["rates"];

	Net_t * net = new Net_t;

	for (size_t i=0; i<inputs.size(); i++)
		net->add_link(inputs[i], rates[i], outputs[i]);

	for (size_t i=0; i<external.size(); i++)
		net->set_source(external[i], ext_rates[i]);

	annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);

	return XPtr<Net_t>(net, true);
	}


size_t find_node_id(const Node_t * n, const Net_t & net)
	{
	for (size_t i=0; i<net.nodes.size(); i++)
		if (n == net.nodes[i])
			return i;

	return ~size_t(0);
	}


void print_PopsNetwork(const XPtr<Net_t> & pNet)
	{
	const Net_t * net = pNet.get();

	Rcout << "Nodes:\n\n";
	Rcout << "id\tinfected\tinput\talleles...\n";
	for (size_t i=0; i<net->nodes.size(); i++)
		{
		const Node_t & n = *net->nodes[i];
		Rcout << i << "\t" <<
			n.rate_in_infd << "\t" <<
			n.rate_in;
		for (auto f : n.frequencies)
			Rcout << "\t" << f;
		Rcout << "\n";
		}
	Rcout << "\n";
	Rcout << "Links:\n\n";
	Rcout << "from\tto\trate\tinfected\n";
	for (size_t i=0; i<net->links.size(); i++)
		{
		L_t & l = *net->links[i];
		size_t f = find_node_id(l.from, *net);
		size_t t = find_node_id(l.to, *net);
		Rcout << f << "\t" <<
			t << "\t" <<
			l.rate << "\t" <<
			l.rate_infd << "\n";
		}
	}


void print_PopsNode(const Node_t * n)
	{
	}

void print_PopsNode(const XPtr<Node_t> & pNode)
	{
	const Node_t * node = pNode.get();
	
	print_PopsNode(node);
	}


XPtr<Net_t> spreadDirichlet(const XPtr<Net_t> & pNet, const List iniDist, double theta)
	{
	Net_t * net = pNet->clone();
	const IntegerVector nodes = iniDist["nodes"];
	const NumericMatrix freqs = iniDist["frequencies"];

	const size_t n_all = freqs.ncol();

	if (nodes.size() != freqs.nrow())
		stop("Invalid parameter 'iniDist': "
		"number of rows in $frequencies and number of elements in $nodes need to be equal!");
	
	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = nodes[i];
		if (n > net->nodes.size())
			stop("Invalid node index in iniDist$nodes!");

		net->nodes[n]->frequencies.resize(n_all);
		for (size_t j=0; j<n_all; j++)
			net->nodes[n]->frequencies[j] = freqs(i, j);
		}

	Drift drift(n_all, theta);
	annotate_frequencies(net->nodes.begin(), net->nodes.end(), drift);
	
	return XPtr<Net_t>(net, true);
	}


XPtr<Node_t> getPopsNode(const XPtr<Net_t> & pNet, int id)
	{
	Net_t * net = pNet.get();
	if (!net)
		stop("Invalid network object!");

	if (size_t(id) > net->nodes.size())
		stop("Invalid node id!");
	
	Node_t * node = net->nodes[size_t(id)];

	// don't GC, since net owns the memory
	return XPtr<Node_t>(node, false);
	}


struct RRng
	{
	double outOf(double mi, double ma) const
		{
		return R::runif(mi, ma);
		}

	size_t operator()(size_t n) const
		{
		size_t r;

		while ((r = R::runif(0, n)) >= n);

		return r;
		}
	};


void sample_node(const Node_t & node, size_t n, vector<size_t> & count)
	{
	if (count.size() != node.frequencies.size())
		stop("Invalid node found!");

	ProportionalPick<> pick(0.000001, node.frequencies);
	RRng r;

	for (size_t i=0; i<n; i++)
		count[pick.pick(r)]++;
	}


IntegerVector drawIsolates_PopsNode(const XPtr<Node_t> & pNode, int n)
	{
	const Node_t * node = pNode.get();
	if (!node)
		stop("Invalid node object!");

	vector<size_t> count(node->frequencies.size(), 0);
	sample_node(*node, n, count);

	return IntegerVector(count.begin(), count.end());
	}


DataFrame drawIsolates_PopsNetwork(const XPtr<Net_t> & pNet, DataFrame samples)
	{
	const Net_t * net = pNet.get();
	if (!net)
		stop("Invalid network object!");

	const IntegerVector nodes = samples["nodes"];
	const IntegerVector num = samples["N"];
	
	// TODO this relies on all nodes having a full frequencies vector!
	const size_t n_freq = net->nodes[0]->frequencies.size();
	vector<IntegerVector> data(n_freq+1);

	for (size_t i=0; i<n_freq+1; i++)
		data[i] = IntegerVector(nodes.size());

	CharacterVector namevec;
	namevec.push_back("node");
	string namestem = "allele_";

	vector<size_t> count(n_freq, 0);

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = nodes[i];
		if (n >= net->nodes.size())
			stop("Invalid node id!");

		sample_node(*net->nodes[n], num[i], count);

		for (size_t j=0; j<n_freq; j++)
			data[j+1][i] = count[j];

		data[0][i] = n;

		namevec.push_back(namestem + to_string(i));
		}

	List dataf(data.begin(), data.end());
	dataf.attr("names") = namevec;
	Rcpp::DataFrame dfout(dataf);
	return dfout;
	}

