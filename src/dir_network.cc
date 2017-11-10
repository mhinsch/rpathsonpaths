#include "dir_network.h"

#include "net_util.h"
#include "rcpp_util.h"
#include "rnet_util.h"
#include "libpathsonpaths/ibmmixed.h"

#include <algorithm>
#include <bitset>


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
	double transmission, double decay, const string & spread_model, bool checks)
	{
	if (checks)
		{
		R_ASSERT(!as<bool>(cycles(links)), "Cycles in network detected");

		IntegerVector subn = colour_network(links);
		R_ASSERT(subn.size() > 0, "Empty network");

		const int col = subn[0];
		for (int c : subn)
			R_ASSERT(c == col, "More than one network in data");
		}

	Net_t * net = new Net_t;

	R_ASSERT(links.size() > 1, 
		"At least two columns required in parameter 'links'.");

	const IntegerVector inputs = links(0);
	const IntegerVector outputs = links(1);
	// this could be done slightly more efficiently but this looks way nicer
	const NumericVector rates = links.size() > 2 ? links(2) :
		NumericVector(inputs.size(), 1.0);

	R_ASSERT(inputs.size() != 0, "Empty network.");
	
	R_ASSERT(external.size() > 1, "At least two columns required in parameter 'external'."); 
	
	const IntegerVector ext_nodes = external(0);
	const NumericVector ext_rates_infd = external(1);
	const NumericVector ext_rates_inp = external.size() > 2 ? external(2) : NumericVector();
	const bool has_inp_rates = ext_rates_inp.size() > 0;

	R_ASSERT(ext_nodes.size() != 0, "No external inputs provided.");

	const int f = int(inputs.inherits("factor")) + outputs.inherits("factor") + 
		ext_nodes.inherits("factor");

	R_ASSERT(f==0 || f==3, "All node lists have to be of the same type.");

	for (auto i : ext_rates_infd)
		R_ASSERT(i >= 0 && i <= 1, "Invalid infected rate.");

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
		try {
		for (size_t i=0; i<ext_nodes.size(); i++)
			if (has_inp_rates)
				net->set_source(ext_nodes[i], ext_rates_infd[i], ext_rates_inp[i]);
			else
				net->set_source(ext_nodes[i], ext_rates_infd[i]);
		} catch (runtime_error & e) {
			stop("Invalid node id in input specification.");
			}

	for (const auto & n : net->nodes)
		R_ASSERT(n != 0, "Invalid network, nodes missing.");

	if (decay >= 0.0 && decay < 1.0)
		preserve_mass(net->nodes.begin(), net->nodes.end(), decay);

	if (spread_model== "fluid")
		// TODO maybe factor out, make constructor only build the net
		annotate_rates(net->nodes.begin(), net->nodes.end(), transmission);
	else if (spread_model == "units")
		{
		Rng rng;
		annotate_rates_ibmm(net->nodes.begin(), net->nodes.end(), transmission, rng);
		}
	else
		stop("Unknown spread model.");

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


XPtr<Net_t> set_allele_freqs(const XPtr<Net_t> & p_net, const List & iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	_set_allele_freqs(net, iniDist);

	return make_S3XPtr(net, "popsnetwork", true);
	}


XPtr<Net_t> popgen_dirichlet(const XPtr<Net_t> & p_net, double theta, Nullable<List> iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	if (! iniDist.isNull())
		_set_allele_freqs(net, iniDist.as());

	R_ASSERT(net->nodes.size(), "Empty network.");

//	print_popsnetwork(make_S3XPtr(net, "popsnetwork", true));

	const size_t n_all = net->nodes[0]->frequencies.size();

	Drift drift(theta);
	annotate_frequencies(net->nodes.begin(), net->nodes.end(), drift);
	
	return make_S3XPtr(net, "popsnetwork", true);
	}


XPtr<Net_t> popgen_ibm_mixed(const XPtr<Net_t> & p_net, Nullable<List> iniDist)
	{
	Net_t * net = new Net_t(*p_net.checked_get());

	if (! iniDist.isNull())
		_set_allele_freqs(net, iniDist.as());

	R_ASSERT(net->nodes.size(), "Empty network");

//	print_popsnetwork(make_S3XPtr(net, "popsnetwork", true));

	const size_t n_all = net->nodes[0]->frequencies.size();

	Rng rng;
	freq_to_popsize_ibmm(net->nodes.begin(), net->nodes.end(), rng);
	annotate_frequencies_ibmm(net->nodes.begin(), net->nodes.end(), rng);
	for (auto node : net->nodes)
		node->normalize();
	
	return make_S3XPtr(net, "popsnetwork", true);
	}


DataFrame draw_isolates_popsnetwork(const XPtr<Net_t> & p_net, const DataFrame & samples)
	{
	const Net_t * net = p_net.checked_get();
	R_ASSERT(net && net->nodes.size()>0, "Invalid or empty network object");

	const IntegerVector nodes = samples(0);
	const IntegerVector num = samples(1);
	const bool f = nodes.inherits("factor");
	const StringVector levels = f ? nodes.attr("levels") : StringVector();

	const size_t n_freq = net->nodes[0]->frequencies.size();
	R_ASSERT(n_freq, "Empty node detected");

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
		R_ASSERT (n < net->nodes.size(), "Invalid node id");

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


DataFrame draw_alleles_popsnetwork
	(const XPtr<Net_t> & p_net, const IntegerVector & nodes, int n)
	{
	const Net_t * net = p_net.checked_get();
	R_ASSERT(net && net->nodes.size()>0, "Invalid or empty network object");

	const bool f = nodes.inherits("factor");
	const StringVector levels = f ? nodes.attr("levels") : StringVector();

	// check nodes
	for (auto n : net->nodes)
		R_ASSERT(n->frequencies.size(), "Empty node detected");

// *** prepare return data

	vector<IntegerVector> data(nodes.size());
	for (auto & v : data)
		v = IntegerVector(n);

// *** generate data

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t nid = f ? net->id_by_name.at(string(levels[nodes[i]-1])):
			nodes[i];

		R_ASSERT(nid < net->nodes.size(), "Invalid node id");

		sample_alleles_node(*net->nodes[nid], data[i]);
		}

// *** construct dataframe and return

	CharacterVector namevec(nodes.size(), "");
	for (size_t i=0; i<nodes.size(); i++)
		if (f) 
			namevec[i] = levels[nodes[i]-1];
		else
			namevec[i] = to_string(i);

	List dataf(nodes.size());
	for (size_t i=0; i<nodes.size(); i++)
		dataf(i) = data[i];
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

		R_ASSERT(l, "Missing link detected");

		const size_t f = net->find_node_id(l->from);
		const size_t t = net->find_node_id(l->to);

		R_ASSERT(f!=n_nodes && t!=n_nodes, "Invalid link");

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


// currently not public
int SNP_distance(int g1, int g2)
	{
	return bitset<sizeof(int)*8>(g1 ^ g2).count();
	}


// currently not public
double distance_SNP(const Node_t & n1, const Node_t & n2)
	{
	double d = 0.0;

	for (int i=0; i<n1.frequencies.size(); i++)
		for (int j=0; j<n2.frequencies.size(); j++)
			d += SNP_distance(i, j) * n1.frequencies[i] * n2.frequencies[j];

	return d;
	}


// currently not public
double SNP_distance_pop(const IntegerVector & p1, const IntegerVector & p2)
	{
	double res = 0.0;

	for (int i1 : p1)
		{
		int sum = 0;
		for (int i2 : p2)
			sum += SNP_distance(i1, i2);
		
		// to avoid overflow we already scale a bit here
		res += sum / p2.size();
		}

	return res/p1.size();
	}


NumericMatrix distances_sample(const XPtr<Net_t> & p_net, int n, bool skip_empty)
	{
	const Net_t * net = p_net.checked_get();

	NumericMatrix res(net->nodes.size(), net->nodes.size());

	vector<vector<size_t>> counts(net->nodes.size());
	for (size_t i=0; i<counts.size(); i++)
		{
		if (skip_empty && net->nodes[i]->rate_in_infd <= 0)
			continue;

		counts[i].resize(net->nodes[i]->frequencies.size(), 0);
		sample_node(*net->nodes[i], n, counts[i]);
		}

	for (int i=0; i<net->nodes.size(); i++)
		for (int j=i; j<net->nodes.size(); j++)
			{
			if (counts[j].empty() || counts[i].empty())
				{
				res(i, j) = res(j, i) = NA_REAL;
				continue;
				}

			for (size_t k=0; k<counts[i].size(); k++)
				res(i, j) += double(abs(int(counts[i][k]) - int(counts[j][k]))); 

			res(i, j) /= 2;
			res(j, i) = res(i, j);
			}

	StringVector cn(net->nodes.size()), rn(net->nodes.size());

	if (net->name_by_id.size())
		{
		// StringVector is clearly missing a constructor here
		cn = net->name_by_id;
		rn = net->name_by_id;
		}
	// we need to name cols and rows even for non-factors, otherwise
	// subscripting won't work (0-based vs. 1-based)
	else
		for (size_t i=0; i<net->nodes.size(); i++)
			cn(i) = rn(i) = to_string(i);

	colnames(res) = cn;
	rownames(res) = rn;

	return res;
	}

NumericMatrix distances_freqdist(const XPtr<Net_t> & p_net, bool skip_empty)
	{
	const Net_t * net = p_net.checked_get();

	NumericMatrix res(net->nodes.size(), net->nodes.size());

	for (int i=0; i<net->nodes.size(); i++)
		for (int j=i; j<net->nodes.size(); j++)
			res(i, j) = res(j, i) = 
				skip_empty && 
					(net->nodes[i]->rate_in_infd <= 0 || net->nodes[j]->rate_in_infd <= 0) ? 
				NA_REAL : distance_freq(*net->nodes[i], *net->nodes[j]);

	if (net->name_by_id.size())
		{
		// StringVector is clearly missing a constructor here
		StringVector cn(net->name_by_id.size()), rn(net->name_by_id.size());
		cn = net->name_by_id;
		rn = net->name_by_id;
		colnames(res) = cn;
		rownames(res) = rn;
		}

	return res;
	}

NumericMatrix distances_EHamming(const XPtr<Net_t> & p_net, bool skip_empty)
	{
	const Net_t * net = p_net.checked_get();

	NumericMatrix res(net->nodes.size(), net->nodes.size());

	for (int i=0; i<net->nodes.size(); i++)
		for (int j=i; j<net->nodes.size(); j++)
			res(i, j) = res(j, i) = 
				skip_empty && 
					(net->nodes[i]->rate_in_infd <= 0 || net->nodes[j]->rate_in_infd <= 0) ? 
				NA_REAL : distance_EHamming(*net->nodes[i], *net->nodes[j]);

	StringVector cn(net->nodes.size()), rn(net->nodes.size());

	if (net->name_by_id.size())
		{
		// StringVector is clearly missing a constructor here
		cn = net->name_by_id;
		rn = net->name_by_id;
		}
	// we need to name cols and rows even for non-factors, otherwise
	// subscripting won't work (0-based vs. 1-based)
	else
		for (size_t i=0; i<net->nodes.size(); i++)
			cn(i) = rn(i) = to_string(i);

	colnames(res) = cn;
	rownames(res) = rn;

	return res;
	}


DataFrame generate_PA(int n_nodes, int n_sources, NumericVector m_dist, float zero_appeal, bool
	compact)
	{
	R_ASSERT(n_sources >= 1, "Number of sources has to be >= 1");
	R_ASSERT(n_nodes >= 1, "Number of nodes has to be >= 1");
	R_ASSERT(zero_appeal > 0, "zero_appeal has to be > 0");

	vector<float> weight(n_nodes + n_sources, 0);

	vector<int> from, to;

	from.reserve(n_nodes);
	to.reserve(n_nodes);
	
	for (size_t i=0; i<n_sources; i++)
		weight[i] = zero_appeal;

	double sum = n_sources * zero_appeal;

	ProportionalPick<> pick(0.000001, m_dist);
	RRng r;

	for (size_t i=0; i<n_nodes; i++)
		{
		// index of current node
		const size_t node = i + n_sources;
		// how many inputs
		const size_t n_inp = pick.pick(r) + 1;

		for (size_t j=0; j<n_inp; j++)
			{
			size_t r_inp = r.outOf(0, sum);
			size_t inp = 0;
			// find previous node in weight list
			while (r_inp > weight[inp])
				r_inp -= weight[inp++];

			R_ASSERT(inp < node, "Invalid node");

			from.push_back(inp);
			to.push_back(node);

			// input node gains a connection
			weight[inp]++;
			// and sum increases accordingly
			sum++;
			}

		// new node has 0 outputs
		weight[node] = zero_appeal;
		sum++;
		}


	if (compact)
		{
		// *** remove isolated nodes === make node indices contiguous
		//     we re-use weight to store how much we have to count the index for a given
		//     node down by

		int reduce = 0;

		// find isolated sources
		for (size_t i=0; i<n_sources; i++)
			{
			// sources that are still at zero_appeal have no outputs => isolated
			if (weight[i] == zero_appeal)
				reduce++;
			else
				weight[i] = reduce;
			}

		// regular nodes can't be isolated, but we still have to change their index 
		fill(weight.begin()+n_sources, weight.end(), reduce);

		for (size_t i=0; i<from.size(); i++)
			{
			from[i] -= weight[from[i]];
			to[i] -= weight[to[i]];
			}
		}

	return DataFrame::create(Named("from") = from, Named("to") = to);
	}
