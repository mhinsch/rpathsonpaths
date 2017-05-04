#include "dir_network.h"

#include "libpathsonpaths/proportionalpick.h"

template<class T>
XPtr<T> make_S3XPtr(T * obj, const char * class_name, bool GC = true)
	{
	XPtr<T> xptr(obj, GC);
	xptr.attr("class") = class_name;
	return xptr;
	}

// TODO factors
IntegerVector sources(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list[0];
	const IntegerVector to = edge_list[1];
	
	vector<bool> is_sink;

	for (auto i : to)
		{
		if (i >= is_sink.size())
			is_sink.resize(i, false);
		is_sink[i] = true;
		}

	set<int> scs;

	for (auto i : from)
		{
		if (i >= is_sink.size() || !is_sink[i])
			scs.insert(i);
		}
	
	IntegerVector res(scs.begin(), scs.end());
	return res;
	}

// TODO factors
IntegerVector colour_network(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list[0];
	const IntegerVector to = edge_list[1];

	// colour of nodes
	vector<int> colour;

	int next_col = 1;

	for (size_t i=0; i<from.size(); i++)
		{
		const size_t f = from[i], t = to[i];

		if (max(f, t) >= colour.size())
			colour.resize(max(f, t), 0);

		if (colour[f] == colour[t])
			{
			// not coloured yet, colour them
			if (colour[f] == 0)
				{
				colour[f] = next_col++;
				colour[t] = colour[f];
				}
			// otherwise they are the same colour which is also fine
			}
		else
			{
			// one of them is not coloured => take the other one's colour
			if (colour[f] == 0)
				{
				colour[f] = colour[t];
				continue;
				}
			if (colour[t] == 0)
				{
				colour[t] = colour[f];
				continue;
				}
			// two different colours, have to change all instances of one of them
			const int oldc = colour[t];
			const int newc = colour[f];

			for (int & c : colour)
				if (c == oldc)
					c = newc;
			}
		}

	IntegerVector res(from.size());

	// assign colours to edges
	for (size_t i=0; i<from.size(); i++)
		res[i] = colour[from[i]];

	return res;
	}


XPtr<Net_t> popsnetwork(const DataFrame & links, const DataFrame & external, 
	double transmission)
	{
	Net_t * net = new Net_t;

	const IntegerVector inputs = links["inputs"];
	const IntegerVector outputs = links["outputs"];
	// this could be done slightly more efficiently but this looks way nicer
	const NumericVector rates = links.containsElementNamed("rates") ? links["rates"] :
		NumericVector(inputs.size(), 1.0);

	const IntegerVector ext_nodes = external["nodes"];
	const NumericVector ext_rates = external["rates"];

	// *** plain integer vector: use ids as raw indices
	if (inputs.attr("class") == "integer")
		{
		const size_t ni = inputs.size();
		for (size_t i=0; i<ni; i++)
			net->add_link(inputs[i], outputs[i], rates[i]);

		if (ext_nodes.attr("class") != "integer")
			stop("Inputs, outputs and external nodes have to be of the same type!");

		const size_t ei = external.size();

		for (size_t i=0; i<ei; i++)
			net->set_source(ext_nodes[i], ext_rates[i]);
		}
	// factor: indices probably differ between input and output, therefore we build
	// our own index
	else if (inputs.attr("class") == "factor")
		{
		const StringVector i_levels = inputs.attr("levels");
		const StringVector o_levels = outputs.attr("levels");
		
		for (size_t i=0; i<inputs.size(); i++)
			{
			// R is 1-based but Rcpp is 0-based
			const auto i_name = i_levels(inputs(i)-1);
			const auto o_name = o_levels(outputs(i)-1);
			
			// gives <iterator, insert succeded>
			const auto i_i = net->lev2idx.emplace(i_name, net->lev2idx.size());
			// insertion succeeded => new name => new node
			if (i_i.second)
				net->idx2lev.push_back(string(i_name));
			// index[name]
			const size_t i_idx = i_i.first->second;
			
			// *** same for output
			const auto o_i = net->lev2idx.emplace(o_name, net->lev2idx.size());
			if (o_i.second)
				net->idx2lev.push_back(string(o_name));

			const size_t o_idx = o_i.first->second;

			net->add_link(i_idx, o_idx, rates[i]);
			}

		if (net->nodes.size() != net->idx2lev.size())
			stop("Something went wrong, mismatch between nodes and levels!");

		if (ext_nodes.attr("class") != "factor")
			stop("Inputs, outputs and external nodes have to be of the same type!");

		const StringVector e_levels = ext_nodes.attr("levels");

		for (size_t i=0; i<ext_nodes.size(); i++)
			{
			// find node with name
			const auto it = net->lev2idx.find(string(e_levels(ext_nodes(i)-1)));
			if (it == net->lev2idx.end())
				stop("Unknown node in 'external'!");
			
			net->set_source(it->second, ext_rates(i));
			}
		}
	else
		stop("List of nodes has to be either integers or factors!");

	for (const auto & n : net->nodes)
		if (n == 0)
			stop("Invalid network, node not set!");

	annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);

	return make_S3XPtr(net, "popsnetwork");
	}

void print_node_id(const Net_t * net, size_t i)
	{
	if (net->idx2lev.size())
		Rcout << net->idx2lev[i];
	else
		Rcout << i;
	}

void print_popsnetwork(const XPtr<Net_t> & p_net)
	{
	const Net_t * net = p_net.checked_get();

	Rcout << "Nodes:\n\n";
	Rcout << "id\tinfected\tinput\talleles...\n";
	for (size_t i=0; i<net->nodes.size(); i++)
		{
		const Node_t & n = *net->nodes[i];
		print_node_id(net, i); Rcout  << "\t" <<
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
		print_node_id(net, f); Rcout  << "\t";
		print_node_id(net, t); Rcout << "\t" <<
			l.rate << "\t" <<
			l.rate_infd << "\n";
		}
	}


void print_popsnode(const Node_t * n)
	{
	}

void print_popsnode(const XPtr<Node_t> & p_node)
	{
	const Node_t * node = p_node.checked_get();
	
	print_popsnode(node);
	}


void _set_allele_freqs(Net_t * net, const List & ini)
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

// TODO factors
XPtr<Net_t> set_allele_freqs(const XPtr<Net_t> & p_net, const List & iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	_set_allele_freqs(net, iniDist);

	return make_S3XPtr(net, "popsnetwork", true);
	}

// TODO factors
XPtr<Net_t> spread_dirichlet(const XPtr<Net_t> & p_net, double theta, Nullable<List> iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	if (! iniDist.isNull())
		_set_allele_freqs(net, iniDist.as());

	if (!net->nodes.size())
		stop("Error: empty network!");

	const size_t n_all = net->nodes[0]->frequencies.size();

	Drift drift(theta);
	annotate_frequencies(net->nodes.begin(), net->nodes.end(), drift);
	
	return make_S3XPtr(net, "popsnetwork", true);
	}

// TODO factor
XPtr<Node_t> get_popsnode(const XPtr<Net_t> & p_net, int id)
	{
	Net_t * net = p_net.checked_get();
	if (!net)
		stop("Invalid network object!");

	if (size_t(id) > net->nodes.size())
		stop("Invalid node id!");
	
	Node_t * node = net->nodes[size_t(id)];

	// don't GC, since net owns the memory
	return make_S3XPtr(node, "popsnode", false);
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


IntegerVector draw_isolates_popsnode(const XPtr<Node_t> & p_node, int n)
	{
	const Node_t * node = p_node.checked_get();
	if (!node)
		stop("Invalid node object!");

	vector<size_t> count(node->frequencies.size(), 0);
	sample_node(*node, n, count);

	return IntegerVector(count.begin(), count.end());
	}


// TODO factors
DataFrame draw_isolates_popsnetwork(const XPtr<Net_t> & p_net, const DataFrame & samples)
	{
	const Net_t * net = p_net.checked_get();
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

// TODO factors
DataFrame edge_list(const XPtr<Net_t> & p_net)
	{
	const Net_t * net = p_net.checked_get();

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

// TODO factors
DataFrame node_list(const XPtr<Net_t> & p_net)
	{
	const Net_t * net = p_net.checked_get();

	StringVector id;
	NumericVector inf;

	for (size_t i=0; i<net->nodes.size(); i++)
		{
		const Node_t * n = net->nodes[i];

		id.push_back(to_string(i));
		inf.push_back(n->rate_in_infd);
		}

	return DataFrame::create(Named("id") = id, Named("infected") = inf);
	}



