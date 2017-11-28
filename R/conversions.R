#' @title epic_to_edge_list
#'
#' @description Convert epicontacts to edgelist.
#'
#' @details Extract an edgelist suitable for popsnetwork from an epicontacts object.
#' Note that this will set all rates to 1.
#'
#' @param x An \code{epicontacts} object.
#'
#' @examples
#' if (require("outbreaks") && require("epicontacts")) {
#'  x <- make_epicontacts(linelist = mers_korea_2015$linelist,
#'                        contacts = mers_korea_2015$contacts)
#'  edgelist <- epic_to_edge_list(x)
#' }
epic_to_edge_list <- function(x) {
    if (!inherits(x, "epicontacts")) {
        stop("x is not an epicontacts object")
    }

    ## The strategy here is to convert whatever was in the edge list to
    ## integers, without gaps, i.e. 0:(N - 1) where N is the number of unique
    ## nodes in the graph

    from_to <- c(x$contacts$from, x$contacts$to)
    from_to_int <- as.integer(factor(from_to)) - 1L
    
    out <- as.data.frame(matrix(from_to_int, ncol = 2))
    rates <- rep(1.0, nrow(out))
    out <- cbind.data.frame(out, rates)
    names(out) <- c("inputs", "outputs", "rates")

    return(out)
}

#' @title rpop_to_dibbler
#' 
#' @description convert an rpopsnetwork object to a dibbler object
#'
#' @details This function converts a popsnetwork object into a dibbler object. 
#' @param network A popsnetwork object.
#' @return A dibbler object.
rpop_to_dibbler <- function(network) {
	if (!requireNamespace("dibbler")){
		stop("This function requires the dibbler package.")
	}
	# get data 
	edge_data <- edge_list(network)
	node_data <- node_list(network)

	dibbler::make_dibbler(edge_data, node_data)
}

