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
#' will be automatically detected and will be processed in order of node id.
#' @param ini_infd Input rate of infected units per source nodes.
#' @param ini_freqs Allele frequencies per source node. This has to be a matrix.
#' @param n Number of simulation runs to perform.
#' @param transmission Rate of infection in nodes.
#' @param decay Decay rate in nodes used to infer transfer rates. If this is negative all 
#' transfer rates will be assumed to be 1 (which in most cases is not a very useful value)
#' unless edgelist has a third column (see above). If decay is not negative transfer rates
#' will be inferred assuming preservation of mass with decay at the nodes. Note that the ibm
#' simulation requires output to be smaller or equal to input for all nodes.
#' @param theta The shape parameter for the Dirichlet distribution.
#' @param method Either "ibm" or "dirichlet".
#' @param checks Whether to perform some sanity checks on the graph before simulating (slow).
#' @return A list containing the result of the simulation(s) as first and the raw network as
#' the second element.
run_popsnet <- function(edgelist, ini_input, ini_infd, ini_freqs, n=1L, transmission=0.0, 
						decay=-1.0, theta=1.0, method="ibm", checks=FALSE) {

	ext_sources <- sort(sources(edgelist))
	n_sources <- length(ext_sources)

	vinp <- rep_len(ini_input, n_sources)
	vinfd <- rep_len(ini_infd, n_sources)

	net_raw <- popsnetwork(edgelist, 
				data.frame(ext_sources, vinfd, vinp), transmission, decay, checks)

	if (n_sources != nrow(ini_freqs)){
		stop(paste("One vector of allele frequencies per input required (got ", 
				   nrow(ini_freqs), " instead of ", n_sources, ")!"))
	}

	if (method == "ibm"){
		res <- replicate(n, spread_ibm_mixed(net_raw, list(ext_sources, ini_freqs))) }
	else if (method == "dirichlet"){
		res <- replicate(n, 
			spread_dirichlet(net_raw, ini_dist=list(ext_sources, ini_freqs), theta=theta)) }
	else {
		stop("Unknown method!") }

	list(result=res, raw=net_raw)
}
