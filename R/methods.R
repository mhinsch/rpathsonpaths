print.PopsNetwork <- function(x, ...){
	.printPopsNetwork(x)
	}

print.PopsNode <- function(x, ...){
	.printPopsNode(x)
	}

drawIsolates <- function(obj, ...){
	UseMethod("drawIsolates")
	}

plot.PopsNetwork <- function(x, ...){
	if (!requireNamespace("igraph")){
		stop("This function require the iGraph package.")
	}

	edges <- edgeList(x)
	graph <- igraph::graph_from_data_frame(edges)
	igraph::plot.igraph(graph)
	}
