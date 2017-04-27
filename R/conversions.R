#' Test file for conversions
#'
#' @param x An \code{epicontacts} object.
#'
#' @export
#'
#' @examples
#' if (require("outbreaks") && require("epicontacts")) {
#'  x <- make_epicontacts(linelist = mers_korea_2015$linelist,
#'                        contacts = mers_korea_2015$contacts)
#'
#' 
#' }
convert <- function(x) {
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
