#'@export
summary.sammon <- function(object,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(object$points, 4))
    cat("\n\n")
    }
