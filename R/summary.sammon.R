#'@export
summary.sammon <- function(object,...)
    {
        out <- object$points
        class(out) <- "summary.sammon"
        out
    }

#' @export
#' @rdname print.sammon
print.summary.sammon <- function(x, digits = max(4L, getOption("digits") - 4L),...)
{
    cat("\n")
    cat("Configurations:\n")
    print(round(x, digits = digits))
    cat("\n")
    invisible(x)
    }
