#' @title print.popsnetwork
#'
#' @description Print a popsnetwork object.
#'
#' @details This will print a textual summary of a popsnetwork object. It will produce two 
#' tables: 
#' \itemize{ 
#' \item A list of all nodes with id (in the same format as provided in the constructor), 
#' proportion
#' of infected material present, overall input rate and genetic composition. 
#' \item A list of all edges with ids of connected nodes, transfer rate and proportion of
#' infected material being transferred.}
#' @param x A popsnetwork object.
print.popsnetwork <- function(x, ...){
	.printpopsnetwork(x)
	}

print.popsnode <- function(x, ...){
	.printpopsnode(x)
	}

draw_isolates <- function(obj, ...){
	UseMethod("draw_isolates")
	}

plot.popsnetwork <- function(x, ...){
	if (!requireNamespace("igraph")){
		stop("This function requires the iGraph package.")
	}

	edges <- edge_list(x)
	colnames(edges)[3] <- "weight"
	max_w <- max(edges[3])
	nodes <- node_list(x)
	graph <- igraph::graph_from_data_frame(edges, TRUE, nodes)
	igraph::E(graph)$width <- igraph::E(graph)$weight/max_w*3
	igraph::plot.igraph(graph, layout=igraph::layout_with_lgl, vertex.size=5, 
						vertex.label.cex=0.5)
	}

#' @title run_popsnet
#'
#' @description Create a popsnetwork and run the transport simulation.
#'
#' @details This function creates a popsnetwork object from the given parameters and runs 
#' either the ibm or dirichlet simulation on it.
#' @param edgelist The list of edges as a dataframe, with origin in the first and target in 
#' the second column. If a third column is present it is interpreted as transfer rates per
#' edge.
#' @param ini_input Net input rates per source node. Note that source nodes in the network
#' will be automatically detected and will be processed in order of node id. ini_input can 
#' either be a single number that will be a repeated for each source node or a list containing
#' one element per source node.
#' @param ini_infd Input rate of infected units per source nodes. Format requirements are
#' identical to those for \code{\link{ini_input}}.
#' @param ini_freqs Allele frequencies per source node. This has to be either a matrix with
#' one row per source node or a list of initialization objects (each of them containing
#' a vector of node ids and a matrix of frequencies, see \code{\link{set_allele_freqs}}). In
#' the latter case n simulations are run per initialization.
#' @param n Number of simulation runs to perform.
#' @param transmission Rate of infection in nodes.
#' @param decay Decay rate in nodes used to infer transfer rates. If this is negative all 
#' transfer rates will be assumed to be 1 (which in most cases is not a very useful value)
#' unless edgelist has a third column (see above). If decay is not negative transfer rates
#' will be inferred assuming preservation of mass with decay at the nodes. Note that the ibm
#' simulation requires output to be smaller or equal to input for all nodes.
#' @param theta The shape parameter for the Dirichlet distribution.
#' @param spread_model Either "fluid" or "units".
#' @param drift_model Either "dirichlet" or "units".
#' @param checks Whether to perform some sanity checks on the graph before simulating (slow).
#' @return A list containing the result of the simulation(s) as first and the raw network as
#' the second element.
run_popsnet <- function(edgelist, ini_input, ini_infd, ini_freqs, n=1L, transmission=0.0, 
						decay=-1.0, theta=1.0, spread_model="units", drift_model="units", 
						checks=FALSE) {

	ext_sources <- sort(sources(edgelist))
	n_sources <- length(ext_sources)

	vinp <- rep_len(ini_input, n_sources)
	vinfd <- rep_len(ini_infd, n_sources)

	net_raw <- popsnetwork(edgelist, 
				data.frame(ext_sources, vinfd, vinp), transmission, decay, checks, 
				spread_model=spread_model)

	if (is.matrix(ini_freqs)){
		if (n_sources != nrow(ini_freqs)){
			stop(paste("One vector of allele frequencies per input required (got ", 
					   nrow(ini_freqs), " instead of ", n_sources, ")!"))
		}

		ini_freqs <- list(list(ext_sources, ini_freqs))
	}

	ini_freqs <- rep(ini_freqs, n)

	if (drift_model== "units"){
		res <- sapply(ini_freqs, function(f) popgen_ibm_mixed(net_raw, f)) }
	else if (drift_model == "dirichlet"){
		res <- sapply(ini_freqs, function(f) popgen_dirichlet(net_raw, theta, f)) }
	else {
		stop("Unknown method!") }

	list(result=res, raw=net_raw)
}

#' @title nodes
#' 
#' @description return a list of nodes
#' 
#' @details Returns a list of node ids for a graph.
#' @param edgelist An edgelist as two columns of node ids (from, to).
#' @return A list or vector of node ids.
nodes <- function(edgelist) unique(c(edgelist[[1]], edgelist[[2]]))

#' @title mutations
#' 
#' @description generate a list of random mutations on a given graph
#' 
#' @details This function will generate a list of mutations suitable as gene frequencies
#' with e.g. \code{\link{set_allele_freqs}} or \code{\link{run_popsnet}}. Each mutation 
#' will be assigned will be assigned to a random node. Allele frequencies will be initialized
#' for the wild type (0) and the mutant (1) and set to 0 otherwise.
#' @param edgelist Two columns of nodes (from, to).
#' @param n_alleles How many alleles to initialize the frequency vector to.
#' @param freq_mutant Frequency of the mutant allele. Note that the wildtype (allele 0)
#' will have frequency 1-freq_mutant.
#' @param n_muts How many mutations to generate.
#' @return Produces a list of mutations. Each mutation consists of a list of node ids and
#' a matrix containing allele frequencies for each node.
mutations <- function(edgelist, n_alleles, freq_mutant, n_muts){
	nods <- nodes(edgelist)
	
	freq <- matrix(rep(0, n_alleles), nrow=1)
	freq[1, 1] <- 1-freq_mutant
	freq[1, 2] <- freq_mutant

	replicate(n_muts, list(sample(nods, 1), freq), simplify=FALSE)
}

descendants <- function(edgelist, node){
	from <- c()
	to <- c()

	use_factor <- is.factor(edgelist[[1]])

	if (is.factor(node)) {
		nodes <- as.character(node)
	} else {
		nodes <- node
	}

	repeat {
		desc <- edgelist[edgelist[1] == nodes,]

		if (nrow(desc) == 0)
			break;

		if (use_factor){
			d1 <- as.character(desc[[1]])
			d2 <- as.character(desc[[2]])
		} else {
			d1 <- desc[[1]]
			d2 <- desc[[2]]
		}

		from <- c(from, d1)
		to <- c(to, d2)

		nodes <- unique(d2)
	}

	if (use_factor){
		from <- factor(from)
		to <- factor(to)
	}

	data.frame(from, to)
}
