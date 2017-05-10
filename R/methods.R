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
	nodes <- node_list(x)
	graph <- igraph::graph_from_data_frame(edges, TRUE, nodes)
	igraph::plot.igraph(graph)
	}
