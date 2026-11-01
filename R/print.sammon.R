#' Print methods for 'sammon' objects
#' 
#' @export
#' @param x an object of class \code{"sammon"} or \code{"summary.sammon"}
#' @param digits integer; number of significant digits used for printing.
#' @param ... further arguments passed to or from other methods.
print.sammon <- function(x, digits = max(2L, getOption("digits") - 2L), ...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: Sammon Scaling \n")
    cat("Number of objects:", dim(x$points)[1], "\n")
    cat("Stress-1:", round(x$stress,digits = digits), "\n")
    cat("\n")
    invisible(x)
    }
