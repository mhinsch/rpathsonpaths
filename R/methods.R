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
		stop("This function require the iGraph package.")
	}

	edges <- edge_list(x)
	nodes <- node_list(x)
	graph <- igraph::graph_from_data_frame(edges, TRUE, nodes)
	igraph::plot.igraph(graph)
	}
