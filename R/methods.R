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
#' @param ... Ignored for now.
#'
#' @examples
#' from <- c(0L, 0L, 1L, 2L, 3L, 1L, 5L)
#' to <- c(1L, 2L, 3L, 3L, 4L, 4L, 2L)
#' rates <- c(1, 1.5, 0.5, 0.1, 1, 0.1, 0.5)
#' edgelist <- data.frame(from, to, rates)
#' ext <- data.frame(c(0L, 5L), c(0.5, 0.5))
#' net <- popsnetwork(edgelist, ext)
#' print(net)
print.popsnetwork <- function(x, ...){
	.printpopsnetwork(x)
	}


#' @title plot.popsnetwork
#'
#' @description Plot a popsnetwork object.
#'
#' @details This function will plot a popsnetwork object as a simple graph (using 
#' \code{\link[igraph]{plot.igraph}}).
#' @param x A popsnetwork object.
#' @param ... Ignored for now.
#' @return NULL
#'
#' @examples
#' from <- c(0L, 0L, 1L, 2L, 3L, 1L, 5L)
#' to <- c(1L, 2L, 3L, 3L, 4L, 4L, 2L)
#' rates <- c(1, 1.5, 0.5, 0.1, 1, 0.1, 0.5)
#' edgelist <- data.frame(from, to, rates)
#' ext <- data.frame(c(0L, 5L), c(0.5, 0.5))
#' net <- popsnetwork(edgelist, ext)
#' plot(net)
plot.popsnetwork <- function(x, ...){
	if (!requireNamespace("igraph")){
		stop("This function requires the iGraph package.")
	}

	edges <- edge_list(x, TRUE)
	colnames(edges)[3] <- "weight"
	max_w <- max(edges[3])
	nodes <- node_list(x, TRUE)
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
#' identical to those for \code{ini_input}.
#' @param ini_freqs Allele frequencies per source node. This can be one out of:
#' * a matrix with one row per root node 
#' * a list of initialization objects (each of them containing a vector of node ids and a matrix of frequencies, see \code{\link{set_allele_freqs}}). In this case n simulations are run per initialization.
#' @param n Number of simulation runs to perform.
#' @param transmission Rate of infection in nodes.
#' @param decay Decay rate in nodes - used to infer transfer rates. If this is negative all 
#' transfer rates will be assumed to be 1 (which in most cases is not a very useful value)
#' unless edgelist has a third column (see above). If decay is not negative transfer rates
#' will be inferred assuming preservation of mass with decay at the nodes. Note that the ibm
#' simulation requires output to be smaller or equal than input for all nodes.
#' @param theta The shape parameter for the Dirichlet distribution.
#' @param spread_model Either "fluid" or "units".
#' @param drift_model Either "dirichlet" or "units".
#' @param checks Whether to perform some sanity checks on the graph before simulating (slow).
#' @return A list containing the network object(s) produced by the simulation(s) as first and 
#' the raw network as the second element.
#'
#' @examples
#' # the network
#' edgelist <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
#' # allele frequencies, 3 alleles x 2 source nodes (A, B)
#' freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
#' # run 
#' run_popsnet(edgelist, 10, 0.1, freqs)
run_popsnet <- function(edgelist, ini_input, ini_infd, ini_freqs, n=1L, transmission=0.0, 
						decay=-1.0, theta=1.0, spread_model="units", drift_model="units", 
						checks=FALSE) {

	meta <- list(spread = spread_model, drift = drift_model, n_inp = length(ini_input))

	ext_sources <- sort(sources(edgelist))
	n_sources <- length(ext_sources)

	meta$n_sources <- n_sources

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

	meta$n_ini <- length(ini_freqs)

	ini_freqs <- rep(ini_freqs, n)

	if (drift_model== "units"){
		res <- sapply(ini_freqs, function(f) popgen_ibm_mixed(net_raw, f)) }
	else if (drift_model == "dirichlet"){
		res <- sapply(ini_freqs, function(f) popgen_dirichlet(net_raw, theta, f)) }
	else {
		stop("Unknown method!") }


	list(result=res, raw=net_raw, meta = meta)
}

#' @title nodes
#' 
#' @description return a list of nodes
#' 
#' @details Returns a list of node ids for a graph.
#' @param edgelist An edgelist as two columns of node ids (from, to).
#' @return A list or vector of node ids.
#'
#' @examples
#' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
#' nodes(el)
nodes <- function(edgelist) {
	if (is.factor(edgelist[[1]])){
		factor(unique(c(as.character(edgelist[[1]]), as.character(edgelist[[2]]))))
	} else {
		unique(c(edgelist[[1]], edgelist[[2]]))
	}
}

#' @title mutations
#' 
#' @description generate a list of random mutations on a given graph
#' 
#' @details This function will generate a list of mutations suitable as gene frequencies
#' with e.g. \code{\link{set_allele_freqs}} or \code{\link{run_popsnet}}. Each mutation 
#' will be assigned to a random node. Allele frequencies will be initialized
#' for the wild type (allele 0) and the mutant (allele 1) and set to 0 otherwise.
#' @param edgelist Two columns of nodes (from, to).
#' @param n_alleles How many alleles to initialize the frequency vector to.
#' @param freq_mutant Frequency of the mutant allele. Note that the wildtype (allele 0)
#' will have frequency \code{1-freq_mutant}.
#' @param n_muts How many mutations to generate.
#' @return Produces a list of mutations. Each mutation consists of a list containing the
#' node id and a 1-row matrix with allele frequencies. This can be used e.g. as the 
#' \code{ini_freqs} parameter of \code{\link{run_popsnet}}.
mutations <- function(edgelist, n_alleles, freq_mutant, n_muts){
	nods <- nodes(edgelist)
	
	freq <- matrix(rep(0, n_alleles), nrow=1)
	freq[1, 1] <- 1-freq_mutant
	freq[1, 2] <- freq_mutant

	replicate(n_muts, list(sample(nods, 1), freq), simplify=FALSE)
}


#' @title descendants
#'
#' @description get all descendants of a particular node in a graph
#' 
#' @details This function will return the graph reachable from a given node in a
#' network. If no node argument is given the function is applied to each node in the
#' network.
#' @param edgelist The graph as an edgelist (a data frame consisting of a from and 
#' a to column).
#' @param node The node to find the descendants of.
#' @return An edgelist containing the complete downstream graph of the given node.
#'
#' @examples
#' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
#' # get descendants of specific nodes
#' descendants(el, "A")
#' descendants(el, c("A", "B"))
#' # apply to all nodes in the network
#' descendants(el)
descendants <- function(edgelist, node){
	if (any(is.na(edgelist))) {
		stop("missing values in edgelist")
	}

	if (missing(node)) {
		ns <- nodes(edgelist)
		return (sapply(ns, function(x) descendants(edgelist, x)))
	}

	from <- c()
	to <- c()

	use_factor <- is.factor(edgelist[[1]])

	if (is.factor(node)) {
		nodes <- as.character(node)
	} else {
		nodes <- node
	}

	repeat {
		desc <- edgelist[edgelist[[1]] %in% nodes,]

		if (nrow(desc) == 0) { break }

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

#' @title biggest_subnetwork
#'
#' @description Find and return the biggest (in terms of number of edges) connected subnetwork
#' of a given network.
#'
#' @details Uses \code{\link{colour_network}} to find isolated subnetworks, then returns only
#' the edges of the subnetwork with the highest number of edges.
#' @param edgelist An edge list in dataframe format.
#' @return An edge list in dataframe format.
biggest_subnetwork <- function(edgelist){
	if (any(is.na(edgelist))) {
		stop("missing values")
	}

	cols <- colour_network(edgelist)
	bigst <- which.max(tabulate(cols))
	edgelist[cols==bigst,]
}

#' @title perfect_binary
#'
#' @description Create a perfect binary tree of a given size.
#'
#' @details Create a perfect, full binary tree with a given depth.
#'
#' @param size Size of the tree in number of edges from the root to any leaf (which is
#' equal to depth-1).
#' @return An edge list in dataframe format.
perfect_binary <- function(size){
	num <- 2 ^ (depth+1) - 2

	from <- vector("integer", num)
	to <- vector("integer", num)

	count <- 1L
	nodes <- c(0L)
	for (i in 1:size){
		# next node id we can use
		next_id <- utils::tail(nodes, 1) + 1L
		# twice as many nodes on this level
		n <- length(nodes) * 2L
		new_nodes <- next_id : (next_id + n-1)
		for (j in 1:n){
			j_o <- (j+1L) %/% 2L
			from[[count]] <- nodes[[j_o]]
			to[[count]] <- new_nodes[[j]]
			count <- count + 1L
		}
		nodes <- new_nodes
	}

	data.frame(from, to)
}

#' @title path_distances
#'
#' @description Obtain the shortest distances between all pairs of nodes.
#'
#' @details Uses \code{\link[igraph]{distances}} to calculate the shortest path for
#' all pairs of nodes in the network.
#'
#' @param net A popsnetwork object.
#' @return The distances in a matrix (see \code{\link[igraph]{distances}}).
path_distances <- function(net) {
	if (!requireNamespace("igraph")){
		stop("This function requires the iGraph package.")
	}

	if (class(net) == "popsnetwork") {
		igraph::distances(igraph::graph_from_data_frame(edge_list(net, TRUE))) 
	} else {
		igraph::distances(igraph::graph_from_data_frame(net)) 
	}
}


#' @title children
#'
#' @description Find all direct children of a node.
#'
#' @details Find and return all children for a given node (or list of nodes) in a network. If 
#' no node is provided the function is applied to all nodes in the network.
#'
#' @param edgelist A network in edgelist format.
#' @param node A single node id or list of node ids.
#' @return An array of nodes or a list of arrays of nodes (if no node was provided).
#'
#' @examples
#' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
#' # get children of specific nodes
#' children(el, "A")
#' children(el, c("A", "B"))
#' # apply to all nodes in the network
#' children(el)
children <- function(edgelist, node) {
	if (any(is.na(edgelist))) {
		stop("missing values in edgelist")
	}

	if (missing(node)) {
		ns <- nodes(edgelist)
		sapply(ns, function(x) children(edgelist, x))
	} else {
		ch <- edgelist[[1]] %in% node
		edgelist[[2]][ch]
	}
}


#' @title parents
#'
#' @description Find all direct parents of a node.
#'
#' @details Find and return all parents for a given node (or list of nodes) in a network. If 
#' no node is provided the function is applied to all nodes in the network.
#'
#' @param edgelist A network in edgelist format.
#' @param node A single node id or list of node ids.
#' @return An array of nodes or a list of arrays of nodes (if no node was provided).
#'
#' @examples
#' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
#' # get parents of specific nodes
#' parents(el, "A")
#' parents(el, c("A", "B"))
#' # apply to all nodes in the network
#' parents(el)
parents <- function(edgelist, node) {
	if (any(is.na(edgelist))) {
		stop("missing values in edgelist")
	}

	if (missing(node)) {
		ns <- nodes(edgelist)
		sapply(ns, function(x) parents(edgelist, x))
	} else {
		ch <- edgelist[[2]] %in% node
		edgelist[[1]][ch]
	}
}


#' @title depth
#' 
#' @description Calculate depth of a node.
#'
#' @details Calculate the depth of a node (or list of nodes) in a given network. The depth is
#' here defined as the length of the shortest path that connects any of the given nodes to
#' any root node. If no node is provided the function is applied to all nodes in the network.
#'
#' @param edgelist The network in edgelist format.
#' @param node A node id or a list of node ids.
#' @return The depth as an integer or an array of integers (if no node was provided).
#'
#' @examples
#' el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"))
#' # get min depth of specific nodes
#' depth(el, "A")
#' depth(el, c("A", "B"))
#' # apply to all nodes in the network
#' # returns a vector, not a number
#' depth(el)
#' # to get the network's depth apply to the list of nodes
#' depth(el, nodes(el))
depth <- function(edgelist, node) {
	if (any(is.na(edgelist))) {
		stop("missing values in edgelist")
	}

	if (missing(node)) {
		ns <- nodes(edgelist)
		sapply(ns, function(x) depth(edgelist, x))
	} else {
		p <- parents(edgelist, node)
		if (length(p) > 0) {
			min(depth(edgelist, p)) + 1L }
		else {
			0L
		}
	}
}
