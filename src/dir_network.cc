#include "dir_network.h"

#include "libpathsonpaths/proportionalpick.h"

template<class T>
XPtr<T> make_S3XPtr(T * obj, const char * class_name, bool GC = true)
	{
	XPtr<T> xptr(obj, GC);
	xptr.attr("class") = class_name;
	return xptr;
	}


XPtr<Net_t> PopsNetwork(const DataFrame & links, const DataFrame & external, double transmission)
	{
	Net_t * net;

	const IntegerVector inputs = links["inputs"];
	const NumericVector rates = links["rates"];
	const IntegerVector outputs = links["outputs"];

	const IntegerVector ext_nodes = external["nodes"];
	const NumericVector ext_rates = external["rates"];

	net = new Net_t;

	const size_t ni = inputs.size();
	const size_t ei = external.size();

	for (size_t i=0; i<ni; i++)
		net->add_link(inputs[i], outputs[i], rates[i]);

	for (size_t i=0; i<ei; i++)
		net->set_source(ext_nodes[i], ext_rates[i]);

	for (const auto & n : net->nodes)
		if (n == 0)
			stop("Invalid network, node not set!");

	annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);

	return make_S3XPtr(net, "PopsNetwork");
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
		Link_t & l = *net->links[i];
		size_t f = net->find_node_id(l.from);
		size_t t = net->find_node_id(l.to);
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


void _setAlleleFreqs(Net_t * net, const List & ini)
	{
	const IntegerVector nodes = ini["nodes"];
	const NumericMatrix freqs = ini["frequencies"];

	const size_t n_all = freqs.ncol();

	if (nodes.size() != freqs.nrow())
		stop("Invalid parameter 'iniDist': "
		"number of rows in $frequencies and number of elements in $nodes need to be equal!");
	
	for (auto n : net->nodes)
		n->frequencies.clear();

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = nodes[i];
		if (n > net->nodes.size())
			stop("Invalid node index in iniDist$nodes!");

		net->nodes[n]->frequencies.resize(n_all, 0);

		for (size_t j=0; j<n_all; j++)
			{
			//Rcout << "setting: " << i << ", " << j << "\n";
			net->nodes[n]->frequencies[j] = freqs(i, j);
			}
		}

	}

XPtr<Net_t> setAlleleFreqs(const XPtr<Net_t> & pNet, const List & iniDist)
	{
	Net_t * net = pNet->clone();

	_setAlleleFreqs(net, iniDist);

	return make_S3XPtr(net, "PopsNetwork", true);
	}

XPtr<Net_t> spreadDirichlet(const XPtr<Net_t> & pNet, double theta, Nullable<List> iniDist)
	{
	Net_t * net = pNet->clone();

	if (! iniDist.isNull())
		_setAlleleFreqs(net, iniDist.as());

	if (!net->nodes.size())
		stop("Error: empty network!");

	const size_t n_all = net->nodes[0]->frequencies.size();

	Drift drift(theta);
	annotate_frequencies(net->nodes.begin(), net->nodes.end(), drift);
	
	return make_S3XPtr(net, "PopsNetwork", true);
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
	return make_S3XPtr(node, "PopsNode", false);
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
		stop("Invalid number of alleles in node!");

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


DataFrame drawIsolates_PopsNetwork(const XPtr<Net_t> & pNet, const DataFrame & samples)
	{
	const Net_t * net = pNet.get();
	if (!net)
		stop("Invalid network object!");

	const IntegerVector nodes = samples["nodes"];
	const IntegerVector num = samples["N"];
	
	const size_t n_freq = net->nodes[0]->frequencies.size();
	if (!n_freq)
		stop("Empty node detected!");

	vector<IntegerVector> data(n_freq+1);

	for (size_t i=0; i<n_freq+1; i++)
		data[i] = IntegerVector(nodes.size());

	CharacterVector namevec;
	namevec.push_back("node");
	string namestem = "allele_";
	for (size_t i=0; i<n_freq; i++)
		namevec.push_back(namestem + to_string(i));

	vector<size_t> count(n_freq, 0);

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = nodes[i];
		if (n >= net->nodes.size())
			stop("Invalid node id!");

		sample_node(*net->nodes[n], num[i], count);

		for (size_t j=0; j<n_freq; j++)
			data[j+1][i] = count[j];

		fill(count.begin(), count.end(), 0);

		data[0][i] = n;
		}

	List dataf(data.begin(), data.end());
	dataf.attr("names") = namevec;
	DataFrame dfout(dataf);
	return dfout;
	}

DataFrame edgeList(const XPtr<Net_t> & pNet)
	{
	const Net_t * net = pNet.get();

	StringVector from;
	StringVector to;
	NumericVector rates;
	NumericVector rates_i;

	const size_t n_nodes = net->nodes.size();

	for (const Link_t * l : net->links)
		{
		if (!l)
			stop("Missing link detected!");

		const size_t f = net->find_node_id(l->from);
		const size_t t = net->find_node_id(l->to);

		if (f==n_nodes || t==n_nodes)
			stop("Invalid link!");

		from.push_back(to_string(f));
		to.push_back(to_string(t));
		rates.push_back(l->rate);
		rates_i.push_back(l->rate_infd);
		}

	return DataFrame::create(
		Named("from") = from,
		Named("to") = to,
		Named("rates") = rates,
		Named("rates_infected") = rates_i);
	}




