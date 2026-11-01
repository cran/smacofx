#' Print methods for 'cmdscalex' objects
#'
#' @export
#' @param x an object of class \code{"cmdscalex"} or \code{"summary.cmdscalex"}
#' @param digits integer; number of significant digits used for printing.
#' @param ... further arguments passed to or from other methods.
print.cmdscalex <- function(x, digits = max(2L, getOption("digits") - 2L), ...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: Torgerson-Gower Scaling \n")
    cat("Number of objects:", dim(x$points)[1], "\n")
    cat("GOF:", round(x$GOF, digits = digits), "\n")
    cat("\n")
    invisible(x)
    }
