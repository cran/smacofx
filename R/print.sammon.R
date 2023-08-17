#'@export
print.sammon <- function(x,...)
    {
    cat("\nCall: ")
    print(x$call)
    cat("\n")
    cat("Model: Sammon Scaling \n")
    cat("Number of objects:", dim(x$points)[1], "\n")
    cat("Stress-1:", round(x$stress,3), "\n")
    cat("\n")
    }
