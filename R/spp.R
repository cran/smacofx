#' Calculating stress per point
#'
#' @param dhat a dist object or symmetric matrix of dissimilarities
#' @param confdist a dist object or symmetric matrix of fitted distances
#' @param weightmat dist objetc or symmetric matrix of weights 
#' 
#' @return a list
#' @export
spp <- function (dhat, confdist, weightmat) 
{
    resmat <- as.matrix(weightmat) * as.matrix(dhat - confdist)^2
    diag(resmat) <- NA
    spp <- colMeans(resmat, na.rm = TRUE)
    spp <- spp/sum(spp) * 100
    names(spp) <- colnames(resmat) <- rownames(resmat) <- attr(dhat, 
        "Labels")
    return(list(spp = spp, resmat = resmat))
}
