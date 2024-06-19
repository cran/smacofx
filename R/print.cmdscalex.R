#'@export
print.cmdscalex <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: Torgerson-Gower Scaling \n")
    cat("Number of objects:", dim(x$points)[1], "\n")
    cat("GOF:", round(x$GOF,2), "\n")
    cat("\n")
    }
