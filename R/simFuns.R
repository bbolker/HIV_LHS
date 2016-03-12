## generic (?) functions for manipulating parameter sets and
## simulation output
## general goals:
##
## * parameter lists:
##     * mutate/transform parameter lists (or should 'mutate' be 'update'?)
##     * transform parameters (dimensionless/dimensioned, constrained/unconstrained spaces, ...)
##     * handle structured parameter lists in a clean way (see ?relist)
##     * reasonable/abbreviated display
##     * pretty-printing
##     * keep associated metadata (parameter descriptions; parent list; purpose/description of parameter set; sources, including citations??)
##
## * define 'flags' (qualitative settings) as well as 'parameters' ?
## * output:
##     * plot (time series, phase plane, ...)

##' change values within a list or vector
##' @param .data a named list
##' @param ... transformation/mutation specifications
##' @importFrom lazyeval all_dots
##' @export
##' @examples
transform.list <- function(.data,...) {
    e <- eval(substitute(list(...)), .data, parent.frame())
    tags <- names(e)
    inx <- match(tags, names(.data))
    matched <- !is.na(inx)
    if (any(matched)) {
        .data[inx[matched]] <- e[matched]
    }
    res <- if (!all(matched)) {
        c(.data, e[!matched])
    } else .data
    ## preserve original class
    class(res) <- class(.data)
    res
}

##' @rdname transform.list
## do we still need transform?
## do we want the ability to mutate within sub-structures?
mutate_.list <- function (.data, ..., .dots)  {
    dots <- lazyeval::all_dots(.dots, ..., all_named = TRUE)
    res <- dplyr:::mutate_impl(.data, dots)
    class(res) <- class(.data)
    return(res)
}


## substitute values into a list
Lsub <- function(L,subs) {
    n <- names(subs)
    for (n in names(subs))
        L[[n]] <- subs[[n]]
    L
}
