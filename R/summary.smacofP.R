

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

#'@export
print.summary.smacofP <- function(x,...)
    {
    cat("\n")
    cat("Configurations:\n")
    print(round(x$conf, 4))
    cat("\n\n")
    cat("Stress per point (in %):\n")
    print(round(x$sppmat[,2], 2))
    cat("\n")
    }

