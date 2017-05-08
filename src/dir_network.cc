#include "dir_network.h"

#include "libpathsonpaths/proportionalpick.h"

#include <algorithm>

template<class T>
XPtr<T> make_S3XPtr(T * obj, const char * class_name, bool GC = true)
	{
	XPtr<T> xptr(obj, GC);
	xptr.attr("class") = class_name;
	return xptr;
	}

// collect a map of factor names to indices from an R factor
vector<size_t> adapt_factor(const IntegerVector & factor, vector<string> & names, 
	unordered_map<string, size_t> & idxs)
	{
	vector<size_t> nodes;

	const StringVector levels = factor.attr("levels");

	for (size_t f : factor)
		{
		// take into account 1-based indexing in R
		const auto name = levels(f-1);

		// try inserting
		const auto ins_it = idxs.emplace(name, names.size());
		// success => name was new => new node
		if (ins_it.second)
			names.push_back(string(name));

		nodes.push_back(ins_it.first->second);
		}

	return nodes;
	}


// wrap an R edge list in order to allow for somewhat sane handling of factors
class EdgeList
	{
	const IntegerVector & _from_raw, & _to_raw;
	vector<size_t> _from, _to;
	vector<string> _names;
	unordered_map<string, size_t> _idxs;

	size_t _f;

public:
	EdgeList(const IntegerVector & from, const IntegerVector & to)
		: _from_raw(from), _to_raw(to)
		{
		_f = int(from.inherits("factor")) + to.inherits("factor");

		if (_f && _f!=2)
			stop("Both node lists have to be of the same type!");

		if (_f)
			{
			_from = adapt_factor(from, _names, _idxs);
			_to = adapt_factor(to, _names, _idxs);
			}
		}

	vector<string> & names() 
		{
		return _names;
		}

	void print()
		{
		for (size_t i=0; i<_from_raw.size(); i++)
			Rcout << _from_raw(i) << " ";
		Rcout << "\n";
		}

	unordered_map<string, size_t> & idxs()
		{
		return _idxs;
		}

	bool factor() const
		{
		return _f;
		}

	size_t from(size_t i) const
		{
		return _f ? _from[i] : _from_raw(i);
		}

	size_t to(size_t i) const
		{
		return _f ? _to[i] : _to_raw(i);
		}

	size_t index(const string & name) const
		{
		return _idxs.at(name);
		}
	const string & name(size_t idx) const
		{
		return _names[idx];
		}
	};

IntegerVector sources(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list[0];
	const IntegerVector to = edge_list[1];

	EdgeList el(from, to);

	vector<bool> is_sink;

	for (size_t i=0; i<to.size(); i++)
		{
		const size_t n = el.to(i);
		if (n >= is_sink.size())
			is_sink.resize(n+1, false);
		is_sink[n] = true;
		}

	set<size_t> scs;

	for (size_t i=0; i<from.size(); i++)
		{
		const size_t n = el.from(i);
		if (n >= is_sink.size() || !is_sink[n])
			scs.insert(n);
		}

	IntegerVector res(scs.size());

	// EdgeList uses 0-base so we might have to rescale here
	size_t i = 0;
	for (size_t s : scs)
		res(i++) = el.factor() ? s+1 : s;

	if (el.factor())
		{
		res.attr("class") = "factor";
		res.attr("levels") = el.names();
		}

	return res;
	}


IntegerVector colour_network(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list[0];
	const IntegerVector to = edge_list[1];

	EdgeList el(from, to);

	// colour of nodes
	vector<int> colour;

	int next_col = 1;

	for (size_t i=0; i<from.size(); i++)
		{
		const size_t f = el.from(i), t = el.to(i);

		if (max(f, t) >= colour.size())
			colour.resize(max(f, t)+1, 0);

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
		res[i] = colour[el.from(i)];

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

	const int f = int(inputs.inherits("factor")) + outputs.inherits("factor") + 
		ext_nodes.inherits("factor");

	if (f!=0 && f!=3)
		stop("All node lists have to be of the same type!");

	EdgeList el(inputs, outputs);

	const size_t ni = inputs.size();
	for (size_t i=0; i<ni; i++)
		net->add_link(el.from(i), el.to(i), rates(i));

	if (el.factor())
		{
		StringVector e_levels = ext_nodes.attr("levels");
		for (size_t i=0; i<ext_nodes.size(); i++)
			net->set_source(el.index(string(e_levels(ext_nodes(i)-1))), ext_rates[i]);

		swap(el.idxs(), net->id_by_name);
		swap(el.names(), net->name_by_id);
		// !!! el is empty below this line !!!
		}
	else
		for (size_t i=0; i<ext_nodes.size(); i++)
			net->set_source(ext_nodes[i], ext_rates[i]);

	for (const auto & n : net->nodes)
		if (n == 0)
			stop("Invalid network, node not set!");

	annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);

	return make_S3XPtr(net, "popsnetwork");
	}

void print_node_id(const Net_t * net, size_t i)
	{
	if (net->name_by_id.size())
		Rcout << net->name_by_id[i];
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


XPtr<Net_t> set_allele_freqs(const XPtr<Net_t> & p_net, const List & iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	_set_allele_freqs(net, iniDist);

	return make_S3XPtr(net, "popsnetwork", true);
	}


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


XPtr<Node_t> get_popsnode(const XPtr<Net_t> & p_net, SEXP id)
	{
	Net_t * net = p_net.checked_get();
	if (!net)
		stop("Invalid network object!");

	const size_t n_id = id_from_SEXP(*net, id);

	if (n_id > net->nodes.size())
		stop("Invalid node id!");
	
	Node_t * node = net->nodes[n_id];

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


DataFrame draw_isolates_popsnetwork(const XPtr<Net_t> & p_net, const DataFrame & samples)
	{
	const Net_t * net = p_net.checked_get();
	if (!net)
		stop("Invalid network object!");

	const IntegerVector nodes = samples["nodes"];
	const IntegerVector num = samples["N"];
	const bool f = nodes.inherits("factor");
	const StringVector levels = f ? nodes.attr("levels") : StringVector();

	const size_t n_freq = net->nodes[0]->frequencies.size();
	if (!n_freq)
		stop("Empty node detected!");

// *** prepare return data

	vector<IntegerVector> data(n_freq+1);
	for (size_t i=0; i<n_freq+1; i++)
		data[i] = IntegerVector(nodes.size());

// *** generate data

	vector<size_t> count(n_freq, 0);

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = f ? net->id_by_name.at(string(levels[nodes[i]-1])):
			nodes[i];
		if (n >= net->nodes.size())
			stop("Invalid node id!");

		sample_node(*net->nodes[n], num[i], count);

		for (size_t j=0; j<n_freq; j++)
			data[j+1][i] = count[j];

		fill(count.begin(), count.end(), 0);

		data[0][i] = n;
		}

// *** construct dataframe and return

	CharacterVector namevec(n_freq+1, "");
	namevec[0] = "node";
	string namestem = "allele_";
	for (size_t i=0; i<n_freq; i++)
		namevec[i+1] = namestem + to_string(i);

	List dataf(n_freq+1);
	dataf(0) = nodes;
	for (size_t i=0; i<n_freq; i++)
		dataf(i+1) = data[i];
	dataf.attr("names") = namevec;

	return DataFrame(dataf);
	}


DataFrame edge_list(const XPtr<Net_t> & p_net)
	{
	const Net_t * net = p_net.checked_get();

	const size_t n_links = net->links.size();

	StringVector from(n_links);
	StringVector to(n_links);
	NumericVector rates(n_links);
	NumericVector rates_i(n_links);

	// do we have names?
	const bool is_factor = net->name_by_id.size();

	const size_t n_nodes = net->nodes.size();

	for (size_t i=0; i<net->links.size(); i++)
		{
		const auto * l = net->links[i];

		if (!l)
			stop("Missing link detected!");

		const size_t f = net->find_node_id(l->from);
		const size_t t = net->find_node_id(l->to);

		if (f==n_nodes || t==n_nodes)
			stop("Invalid link!");

		from[i] = is_factor ? net->name_by_id[f] : to_string(f);
		to[i] = is_factor ? net->name_by_id[t] : to_string(t);
		rates[i] = l->rate;
		rates_i[i] = l->rate_infd;
		}

	return DataFrame::create(
		Named("from") = from,
		Named("to") = to,
		Named("rates") = rates,
		Named("rates_infected") = rates_i);
	}


DataFrame node_list(const XPtr<Net_t> & p_net)
	{
	const Net_t * net = p_net.checked_get();

	StringVector id(net->nodes.size());
	NumericVector inf(net->nodes.size());

	const bool is_factor = net->name_by_id.size();

	for (size_t i=0; i<net->nodes.size(); i++)
		{
		const Node_t * n = net->nodes[i];

		id[i] = is_factor ? net->name_by_id[i] : to_string(i);
		inf[i] = n->rate_in_infd;
		}

	return DataFrame::create(Named("id") = id, Named("infected") = inf);
	}



