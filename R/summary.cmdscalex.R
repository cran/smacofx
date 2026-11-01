#'@export
summary.cmdscalex <- function(object,...)
    {
        out <- object$points
        class(out) <- "summary.cmdscalex"
        out
    }


#' @export
#' @rdname print.cmdscalex
print.summary.cmdscalex <- function(x, digits = max(4L, getOption("digits") - 4L), ...)
{
    cat("\n")
    cat("Configurations:\n")
    print(round(x, digits = digits))
    cat("\n")
    invisible(x)
    }
