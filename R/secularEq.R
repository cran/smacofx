

#' Secular Equation 
#'
#' @param a matrix
#' @param b matrix
#'
#' @importFrom stats uniroot
secularEq<-function(a,b) {
    n<-dim(a)[1]
    eig<-eigen(a)
    eva<-eig$values
    eve<-eig$vectors
    beta<-drop(crossprod(eve, b))
    f<-function(mu) {
        return(sum((beta/(eva+mu))^2)-1)
    }
    lmn<-eva [n]
    uup<-sqrt(sum(b^2))-lmn
    ulw<-abs(beta [n])-lmn
    rot<-stats::uniroot(f,lower= ulw,upper= uup)$root
    cve<-beta/(eva+rot)
    return(drop(eve%*%cve))
}    
