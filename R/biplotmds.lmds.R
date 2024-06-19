#' S3 method for lmds objects
#'
#' @param object An object of class smacofP
#' @param extvar Data frame with external variables.
#' @param scale if 'TRUE' external variables are standardized internally.
#'
#' @details
#' If a model for individual differences is provided, the external
#' variables are regressed on the group stimulus space
#' configurations. For objects returned from 'biplotmds' we use the plot method in
#' \code{\link[smacof]{biplotmds}}. In the biplot called with plot() only the relative length of the
#' vectors and their direction matters. Using the vecscale argument in plot() the
#' user can control for the relative length of the vectors. If
#' 'vecscale = NULL', the 'vecscale()' function from the 'candisc'
#' package is used which tries to automatically calculate the scale
#' factor so that the vectors approximately fill the same space as
#' the configuration. In this method vecscale should usually be smaller than the one used in smacof
#' by a factor of 0.1. Note that in the biplot object, the configuration is always normalized (which it may not necessarily be in the lmds object).
#'
#' @return
#' Returns an object belonging to classes 'mlm' and 'mdsbi'. See 'lm' for details.
#' R2vec: Vector containing the R2 values.
#' See also \code{\link[smacof]{biplotmds}} for the plot method. 
#'
#'
#' @importFrom smacof biplotmds
#'
#'
#' @export
#' 
biplotmds.lmds <- function(object,extvar,scale=TRUE)
{
    oldclass <- class(object)
    if(any(class(object)%in% "smacof")) class(object) <- "smacof"
    oldconf <- object$conf
    object$conf <- object$conf/enorm(object$conf)
    out <- smacof::biplotmds(object,extvar=extvar,scale=scale)
    out$conf <- object$conf
    object$conf <- oldconf
    class(object) <- oldclass
    return(out)
}



