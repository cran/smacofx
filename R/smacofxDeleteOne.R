#' Helper function to conduct jackknife MDS
#'
#' Function deletes every object row and columns once and fits the MDS in object and returns the configuration. The deleted row is set to 0 in the configuration. Is meant for smacofx functions, but should also work for every smacof models. 
#' 
#' @param object  Object of class smacofP if used as method or another object inheriting from smacofB. Note we assume the MDS model was fitted on a symmetric matrix/data frame or dist object 
#' @param delta the data (symmetric matrix, data frame or dist object)
#' @param ndim target dimension of the mds
#' @param type type of MDS
#' @param init starting configuration 
#' @param weightmat weighting matrix
#' @param verbose print progress
#' @param itmaxi maximum iterations of the MDS procedure
#'
#' @return An array of size n with n coonfigurations
smacofxDeleteOne <- function (object, delta, weightmat, init, ndim, type, verbose=FALSE, itmaxi=10000) {
    ocall<-as.call(object$call)
    if(is.numeric(ocall$itmax)) itmaxi <- ocall$itmax #if the call has an itmax agrumnet we use that otherwise we use what was supplied in jackmds.smacofP
    n <- nrow (delta)
    x <- array (0, c (n, ndim, n))
    reslist <- list()
    flag <- any(class(object)%in%c("bcmds","lmds")) 
    for (i in 1:n) {
        if(isTRUE(verbose)) cat("Jackknife Sample: ", formatC(i, digits = 3, width = 3))
        #xi <- smacofSym(delta[-i, -i], ndim = ndim, type = type)$conf
        ocall$delta<-delta[-i, -i]
        if(!flag) ocall$weightmat<-weightmat[-i, -i]
        ocall$init <- init[-i,]
        ocall$itmax <- itmaxi
        xi <- eval(ocall)
        xconf <- xi$conf
        x[((1 : n)[-i]), (1 : ndim), i] <- xconf
        x[i, (1 : ndim), i] <- 0
        reslist[[i]] <- xi
        if(isTRUE(verbose)) cat("\n")      
        }
    return (list(x=x,res=reslist))
}
