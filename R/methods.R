print.PopsNetwork <- function(x, ...){
	.printPopsNetwork(x)
	}

print.PopsNode <- function(x, ...){
	.printPopsNode(x)
	}

drawIsolates <- function(obj, ...){
	UseMethod("drawIsolates")
	}
