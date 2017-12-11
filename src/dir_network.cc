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

	// convert to EdgeList so we don't have to worry about factor/integer
	EdgeList el(from, to);

	const set<size_t> scs = find_sources(el);

	IntegerVector res(scs.size());

	size_t i = 0;
	for (size_t s : scs)
		// EdgeList uses 0-base so we might have to rescale here
		res(i++) = el.factor() ? s+1 : s;

	// build a factor if necessary
	el.make_factor(res);

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

	el.make_factor(res);

	return res;
	}


IntegerVector colour_network(const DataFrame & edge_list)
	{
	const IntegerVector from = edge_list(0);
	const IntegerVector to = edge_list(1);

	EdgeList el(from, to);

	// colour of nodes
	vector<int> colour = colour_network(el.begin(), el.end());

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
		// find cycles for each source node
		for (size_t i : scs)
			cycles.find_cycles(i);
		
		// build an R compatible list (of int vectors or factors) from the result
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
	// do some slow sanity checks
	if (checks)
		{
		// cycles
		R_ASSERT(!as<bool>(cycles(links)), "Cycles in network detected");

		// non-empty
		IntegerVector subn = colour_network(links);
		R_ASSERT(subn.size() > 0, "Empty network");

		// no separate networks
		const int col = subn[0];
		for (int c : subn)
			R_ASSERT(c == col, "More than one network in data");
		}

	Net_t * net = new Net_t;

	R_ASSERT(links.size() > 1, 
		"At least two columns required in parameter 'links'.");

// *** edge list, rates
	const IntegerVector inputs = links(0);		// needed
	const IntegerVector outputs = links(1);		// needed
	// transfer rates are optional (default value is 1)
	// this could be done slightly more efficiently but this looks way nicer
	const NumericVector rates = links.size() > 2 ? links(2) :
		NumericVector(inputs.size(), 1.0);

	R_ASSERT(inputs.size() != 0, "Empty network.");
	R_ASSERT(external.size() > 1, "At least two columns required in parameter 'external'."); 

// *** external input
	const IntegerVector ext_nodes = external(0);		// nodes are needed
	const NumericVector ext_rates_infd = external(1);	// rate of infectedness is needed
	// input rates are optional
	const NumericVector ext_rates_inp = external.size() > 2 ? external(2) : 
		NumericVector(ext_nodes.size(), 1.0);
	const bool has_inp_rates = external.size() > 2;

	R_ASSERT(ext_nodes.size() != 0, "No external inputs provided.");

	const int f = int(inputs.inherits("factor")) + outputs.inherits("factor") + 
		ext_nodes.inherits("factor");
	R_ASSERT(f==0 || f==3, "All node lists have to be of the same type.");

	EdgeList el(inputs, outputs);

// *** basic topology
	const size_t ni = inputs.size();
	for (size_t i=0; i<ni; i++)
		net->add_link(el.from(i), el.to(i), rates(i));

// *** external inputs
	if (el.factor())
		{
		StringVector e_levels = ext_nodes.attr("levels");
		for (size_t i=0; i<ext_nodes.size(); i++)
			{
			R_ASSERT(ext_rates_infd[i] <= ext_rates_inp[i], 
				"input of infected material larger than overall input");
			net->set_source(el.index(string(e_levels(ext_nodes(i)-1))), 
				ext_rates_infd[i], ext_rates_inp[i]);
			}

		// TODO not pretty, should be done better
		swap(el.idxs(), net->id_by_name);
		swap(el.names(), net->name_by_id);
		// !!! el is empty below this line !!!
		}
	else
		try {
		for (size_t i=0; i<ext_nodes.size(); i++)
			{
			R_ASSERT(ext_rates_infd[i] <= ext_rates_inp[i], 
				"input of infected material larger than overall input");
			net->set_source(ext_nodes[i], ext_rates_infd[i], ext_rates_inp[i]);
			}
		} catch (runtime_error & e) {
			stop("Invalid node id in input specification.");
			}

	// check for gaps in ids
	for (const auto & n : net->nodes)
		R_ASSERT(n != 0, "Invalid network, nodes missing.");

	// this interpolates transfer rates
	if (decay >= 0.0 && decay < 1.0)
		preserve_mass(net->nodes.begin(), net->nodes.end(), decay);

// *** generate rate of infectedness for all nodes
	if (spread_model== "fluid")
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

	const size_t n_all = net->nodes[0]->frequencies.size();

	R_ASSERT(n_all, "No genetic data in network.");

	// simulate
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

	const size_t n_all = net->nodes[0]->frequencies.size();

	R_ASSERT(n_all, "No genetic data in network.");

	Rng rng;
	// scale frequencies to absolute numbers
	freq_to_popsize_ibmm(net->nodes.begin(), net->nodes.end(), rng);
	// simulate
	annotate_frequencies_ibmm(net->nodes.begin(), net->nodes.end(), rng);
	// scale back to frequencies
	for (auto node : net->nodes)
		node->normalize();
	
	return make_S3XPtr(net, "popsnetwork", true);
	}


DataFrame draw_isolates(const XPtr<Net_t> & p_net, const DataFrame & samples, bool aggregate)
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

	// if we don't aggregate we have one line per sample
	// otherwise it's one line per node
	const size_t n_data = aggregate ? 
		num.size() : std::accumulate(num.begin(), num.end(), 0);

	vector<IntegerVector> data(aggregate ? n_freq : 1);
	for (auto & i : data) 
		i = IntegerVector(n_data);

	// if we don't aggregate we have to expand the node list by #samples per node
	IntegerVector nodes_r(aggregate ? 0 : n_data);
	if (!aggregate)
		{
		size_t i = 0;
		for (size_t n=0; n<nodes.size(); n++)
			{
			const int node = nodes[n];
			const int n_samples = num[n];
			R_ASSERT(i+n_samples <= n_data, "something seriously wrong");
			for (size_t j=0; j<n_samples; j++)
				nodes_r[i++] = node;
			}
		// convert to factor if necessary
		if (f)
			{
			nodes_r.attr("class") = "factor";
			nodes_r.attr("levels") = levels;
			}
		}

// *** generate data

	vector<size_t> count(n_freq, 0);

	// counter for non-aggregate case
	size_t na_idx= 0;
	vector<int> na_count;

	for (size_t i=0; i<nodes.size(); i++)
		{
		const size_t n = f ? net->id_by_name.at(string(levels[nodes[i]-1])):
			nodes[i];
		R_ASSERT (n < net->nodes.size(), "Invalid node id");


		if (aggregate)
			{
			// draw samples
			sample_node(*net->nodes[n], num[i], count);
			// copy to data frame
			for (size_t j=0; j<n_freq; j++)
				data[j][i] = count[j];
			// reset for next round
			fill(count.begin(), count.end(), 0);
			}
		else
			{
			// size of vector determines #samples
			na_count.resize(num[i], 0);
			// draw samples
			sample_alleles_node(*net->nodes[n], na_count);
			// copy to return vector
			for (int allele : na_count)
				data[0][na_idx++] = allele;
			}
		}

// *** construct dataframe and return

	const size_t n_cols = aggregate ? n_freq+1 : 2;

	CharacterVector namevec(n_cols, "");
	namevec[0] = "node";
	if (aggregate)
		{
		const string namestem = "allele_";
		for (size_t i=0; i<n_freq; i++)
			namevec[i+1] = namestem + to_string(i);
		}
	else
		namevec[1] = "allele";

	List dataf(n_cols);
	if (aggregate)
		{
		dataf(0) = nodes;
		for (size_t i=0; i<n_freq; i++)
			dataf(i+1) = data[i];
		}
	else
		{
		dataf(0) = nodes_r;
		dataf(1) = data[0];
		}

	dataf.attr("names") = namevec;

	return DataFrame(dataf);
	}


DataFrame draw_alleles(const XPtr<Net_t> & p_net, const IntegerVector & nodes, int n)
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


DataFrame edge_list(const XPtr<Net_t> & p_net, bool as_string)
	{
	const Net_t * net = p_net.checked_get();

	const size_t n_links = net->links.size();

	StringVector from_s(as_string ? n_links : 0);
	StringVector to_s(as_string ? n_links : 0);

	IntegerVector from_i(as_string ? 0 : n_links);
	IntegerVector to_i(as_string ? 0 : n_links);

	NumericVector rates(n_links);
	NumericVector rates_i(n_links);

	// do we have names?
	const bool is_factor = net->name_by_id.size();

	const size_t n_nodes = net->nodes.size();

	// fill lists
	for (size_t i=0; i<net->links.size(); i++)
		{
		const auto * l = net->links[i];

		R_ASSERT(l, "Missing link detected");

		const size_t f = net->find_node_id(l->from);
		const size_t t = net->find_node_id(l->to);

		R_ASSERT(f!=n_nodes && t!=n_nodes, "Invalid link");

		if (as_string)
			{
			from_s[i] = is_factor ? net->name_by_id[f] : to_string(f);
			to_s[i] = is_factor ? net->name_by_id[t] : to_string(t);
			}
		else
			{
			from_i[i] = is_factor ? f+1 : f;
			to_i[i] = is_factor ? t+1 : t;
			}

		rates[i] = l->rate;
		rates_i[i] = l->rate_infd;
		}

	if (as_string)
		return DataFrame::create(
			Named("from") = from_s,
			Named("to") = to_s,
			Named("rates") = rates,
			Named("rates_infected") = rates_i,
			// this is a ridiculous R-ism
			Named("stringsAsFactors") = false);

	if (is_factor)
		{
		from_i.attr("class") = "factor";
		from_i.attr("levels") = net->name_by_id;
		to_i.attr("class") = "factor";
		to_i.attr("levels") = net->name_by_id;
		}

	return DataFrame::create(
		Named("from") = from_i,
		Named("to") = to_i,
		Named("rates") = rates,
		Named("rates_infected") = rates_i);
	}


DataFrame node_list(const XPtr<Net_t> & p_net, bool as_string)
	{
	const Net_t * net = p_net.checked_get();

	StringVector id_s(as_string ? net->nodes.size() : 0);
	IntegerVector id_i(as_string ? 0 : net->nodes.size());
	NumericVector inf(net->nodes.size());

	const bool is_factor = net->name_by_id.size();

	for (size_t i=0; i<net->nodes.size(); i++)
		{
		const Node_t * n = net->nodes[i];

		if (as_string)
			id_s[i] = is_factor ? net->name_by_id[i] : to_string(i);
		else
			id_i[i] = is_factor ? i+1 : i;

		inf[i] = n->rate_in_infd;
		}

	if (as_string)
		return DataFrame::create(Named("id") = id_s, Named("infected") = inf);

	if (is_factor)
		{
		id_i.attr("class") = "factor";
		id_i.attr("levels") = net->name_by_id;
		}

	return DataFrame::create(Named("id") = id_i, Named("infected") = inf);
	}


NumericMatrix distances_sample(const XPtr<Net_t> & p_net, int n, bool skip_empty)
	{
	const Net_t * net = p_net.checked_get();

	NumericMatrix res(net->nodes.size(), net->nodes.size());

	R_ASSERT(net->nodes.size(), "empty network detected");

	// allele counts per node
	vector<vector<size_t>> counts(net->nodes.size());
	for (size_t i=0; i<counts.size(); i++)
		{
		if (skip_empty && net->nodes[i]->rate_in_infd <= 0)
			continue;

		R_ASSERT(net->nodes[i]->frequencies.size(), "no genetic data in network");

		counts[i].resize(net->nodes[i]->frequencies.size(), 0);
		// draw n samples from node, count occurence of each allele
		sample_node(*net->nodes[i], n, counts[i]);
		}

	// calculate distances
	for (int i=0; i<net->nodes.size(); i++)
		for (int j=i; j<net->nodes.size(); j++)
			{
			// one or both empty => NA
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

	// col/row names
	StringVector cn(net->nodes.size()), rn(net->nodes.size());

	if (net->name_by_id.size())		// factor
		{
		// StringVector is clearly missing a constructor here
		cn = net->name_by_id;
		rn = net->name_by_id;
		}
	// we need to name cols and rows even for non-factors, otherwise
	// subscripting won't work (0-based vs. 1-based)
	else							// not factor
		for (size_t i=0; i<net->nodes.size(); i++)
			cn(i) = rn(i) = to_string(i);

	colnames(res) = cn;
	rownames(res) = rn;

	return res;
	}

NumericMatrix distances_freqdist(const XPtr<Net_t> & p_net, bool skip_empty)
	{
	const Net_t * net = p_net.checked_get();

	R_ASSERT(net->nodes.size(), "empty network detected");

	NumericMatrix res(net->nodes.size(), net->nodes.size());

	// distances
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

	R_ASSERT(net->nodes.size(), "empty network detected");

	NumericMatrix res(net->nodes.size(), net->nodes.size());

	// distances
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

	vector<int> from, to;
	
	// distribution for #input links
	ProportionalPick<> pick(0.000001, m_dist);
	RRng r;

	net_gen_prefattach(from, to, n_nodes, n_sources, 
		[&pick,&r] (int) -> int {return pick.pick(r);}, 
		zero_appeal, r, compact);

	return DataFrame::create(Named("from") = from, Named("to") = to);
	}
