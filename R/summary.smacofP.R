#'@export
summary.smacofP <- function(object,...)
    {
      spp.perc <- object$spp/sum(object$spp) * 100
      sppmat <- cbind(sort(object$spp), sort(spp.perc))
      colnames(sppmat) <- c("SPP", "SPP(%)") 
      res <- list(conf=object$conf,sppmat=sppmat)
      class(res) <- "summary.smacofP"
      res
    }

#' @export
#' @rdname print.smacofP
print.summary.smacofP <- function(x, digits = max(4L, getOption("digits") - 4L), ...)
{
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, digits = digits))
    cat("\n\n")
    cat("Stress per point (in %):\n")
    print(round(x$sppmat[,2], digits = digits))
    cat("\n")
    invisible(x)
    }

