#include "dir_network.h"

#include "net_util.h"
#include "rcpp_util.h"
#include "rnet_util.h"
#include "libpathsonpaths/ibmmixed.h"

#include <algorithm>


IntegerVector sources(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list(0);
	const IntegerVector to = edge_list(1);

	EdgeList el(from, to);

	const set<size_t> scs = find_sources(el);

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

IntegerVector sinks(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list(0);
	const IntegerVector to = edge_list(1);

	EdgeList el(from, to);

	const set<size_t> scs = find_sinks(el);

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
	const IntegerVector from = edge_list(0);
	const IntegerVector to = edge_list(1);

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


SEXP cycles(const DataFrame & edge_list, bool record)
	{
	const IntegerVector from = edge_list(0);
	const IntegerVector to = edge_list(1);

	EdgeList el(from, to);

	const set<size_t> scs = find_sources(el);
	const size_t n_nodes = el.n_nodes();

	// convert edge list to table
	vector<vector<size_t> > outputs(n_nodes);
	for (size_t i=0; i<from.size(); i++)
		outputs[el.from(i)].push_back(el.to(i));

	Cycles cycles(outputs);

	// user wants a list of cycles
	if (record)
		{
		for (size_t i : scs)
			cycles.find_cycles(i);
		
		List res(cycles.res.size());
		if (el.factor())
			{
			for (size_t i=0; i<cycles.res.size(); i++)
				{
				IntegerVector v(cycles.res[i].size());
				// convert to 1-based indexing
				for (size_t j=0; j<cycles.res[i].size(); j++)
					v[j] = cycles.res[i][j] + 1;
				v.attr("class") = "factor";
				v.attr("levels") = el.names();
				res[i] = v;
				}
			}
		else
			for (size_t i=0; i<cycles.res.size(); i++)
				res[i] = cycles.res[i];

		return res;
		}
	// yes or no is fine
	else
		{
		for (size_t i : scs)
			if (cycles.has_cycles(i))
				return wrap(true);

		return wrap(false);
		}
	}


XPtr<Net_t> popsnetwork(const DataFrame & links, const DataFrame & external, 
	double transmission, bool checks)
	{
	if (checks)
		{
		if (as<bool>(cycles(links)))
			stop("Cycles in network detected!");

		IntegerVector subn = colour_network(links);
		if (subn.size() == 0)
			stop("Empty network!");
		const int col = subn[0];

		for (int c : subn)
			if (c != col)
				stop("More than one network in data!");
		}

	Net_t * net = new Net_t;

	const IntegerVector inputs = links(0);
	const IntegerVector outputs = links(1);
	// this could be done slightly more efficiently but this looks way nicer
	const NumericVector rates = links.size() > 2 ? links(2) :
		NumericVector(inputs.size(), 1.0);

	if (inputs.size() == 0)
		stop("Empty network!");

	const IntegerVector ext_nodes = external(0);
	const NumericVector ext_rates_infd = external(1);
	const NumericVector ext_rates_inp = external.size() > 2 ? external(3) : NumericVector();
	const bool has_inp_rates = ext_rates_inp.size() > 0;

	if (ext_nodes.size() == 0)
		stop("No external inputs provided!");

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
			if (has_inp_rates)
				net->set_source(el.index(string(e_levels(ext_nodes(i)-1))), 
					ext_rates_infd[i], ext_rates_inp[i]);
			else
				net->set_source(el.index(string(e_levels(ext_nodes(i)-1))), 
					ext_rates_infd[i]);

		// TODO not pretty, should be done better
		swap(el.idxs(), net->id_by_name);
		swap(el.names(), net->name_by_id);
		// !!! el is empty below this line !!!
		}
	else
		for (size_t i=0; i<ext_nodes.size(); i++)
			if (has_inp_rates)
				net->set_source(ext_nodes[i], ext_rates_infd[i], ext_rates_inp[i]);
			else
				net->set_source(ext_nodes[i], ext_rates_infd[i]);

	for (const auto & n : net->nodes)
		if (n == 0)
			stop("Invalid network, node not set!");

	// TODO maybe factor out, make constructor only build the net
	annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);

	return make_S3XPtr(net, "popsnetwork");
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
			(n.rate_in <= 0 ? 0 : n.rate_in_infd/n.rate_in) << "\t" <<
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
			(l.rate <= 0 ? 0 : l.rate_infd/l.rate) << "\n";
		}
	}


void print_popsnode(const XPtr<Node_t> & p_node)
	{
	const Node_t * node = p_node.checked_get();
	
	print_popsnode(node);
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


XPtr<Net_t> spread_ibm_mixed(const XPtr<Net_t> & p_net, Nullable<List> iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	if (! iniDist.isNull())
		_set_allele_freqs(net, iniDist.as());

	if (!net->nodes.size())
		stop("Error: empty network!");

	const size_t n_all = net->nodes[0]->frequencies.size();

	Rng rng;
	annotate_frequencies_ibmm(net->nodes.begin(), net->nodes.end(), rng);
	
	return make_S3XPtr(net, "popsnetwork", true);
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


IntegerVector draw_isolates_popsnode(const XPtr<Node_t> & p_node, int n)
	{
	const Node_t * node = p_node.checked_get();
	if (!node)
		stop("Invalid node object!");

	vector<size_t> count(node->frequencies.size(), 0);
	if (!count.size())
		stop("Empty node!");

	sample_node(*node, n, count);

	return IntegerVector(count.begin(), count.end());
	}


DataFrame draw_isolates_popsnetwork(const XPtr<Net_t> & p_net, const DataFrame & samples)
	{
	const Net_t * net = p_net.checked_get();
	if (!net || net->nodes.size()==0)
		stop("Invalid or empty network object!");

	const IntegerVector nodes = samples(0);
	const IntegerVector num = samples(1);
	const bool f = nodes.inherits("factor");
	const StringVector levels = f ? nodes.attr("levels") : StringVector();

	const size_t n_freq = net->nodes[0]->frequencies.size();
	if (!n_freq)
		stop("Empty node detected!");

// *** prepare return data

	vector<IntegerVector> data(n_freq);
	for (size_t i=0; i<n_freq; i++)
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
			data[j][i] = count[j];

		fill(count.begin(), count.end(), 0);
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



