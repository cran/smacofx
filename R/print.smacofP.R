#' Print methods for 'smacofP' objects
#'
#' @export
#' @param x an object of class \code{"smacofP"} or \code{"summary.smacofP"}.
#' @param digits integer; number of significant digits used for printing.
#' @param ... further arguments passed to or from other methods.
print.smacofP <- function(x, digits = max(2L, getOption("digits") - 2L), ...)
{
    cat("\nCall:\n")
    print(x$call)
    cat("\n")
    cat("Model:", x$type, x$model, "with parameter vector=",unlist(x$parameters),"\n")
    cat("Number of objects:", x$nobj, "\n")
    cat("Stress-1 value:", round(x$stress, digits = digits), "\n")
    cat("Number of iterations:", x$niter, "\n")
    cat("\n")
    invisible(x)
}
