#' Box-Cox MDS
#'
#' This function minimizes the Box-Cox Stress of Chen & Buja (2013) via gradient descent. This is a ratio metric scaling method. The transformations are not straightforward to interpret but mu is associated with fitted distances in the configuration and lambda with the dissimilarities. Concretely for fitted distances (attraction part) it is \eqn{BC_{mu+lambda}(d(X))} and for the repulsion part it is \eqn{delta^lambda BC_{mu}(d(X))} with BC being the one-parameter Box-Cox transformation.
#'
#' 
#' @param delta dissimilarity or distance matrix, dissimilarity or distance data frame or 'dist' object
#' @param mu mu parameter. Should be 0 or larger for everything working ok. If mu<0 it works but I find the MDS model is strange and normalized stress tends towards 0 regardless of fit. Use normalized stress at your own risk in that case.
#' @param lambda lambda parameter. Must be larger than 0.
#' @param rho the rho parameter, power for the weights (called nu in the original article).
#' @param type what type of MDS to fit. Only "ratio" currently. 
#' @param ndim the dimension of the configuration
#' @param init initial configuration. If NULL a classical scaling solution is used. 
#' @param weightmat a matrix of finite weights. Not implemented.
#' @param itmax number of optimizing iterations, defaults to 2000.
#' @param verbose prints progress if > 3.
#' @param acc Accuracy (lowest stepsize). Defaults to 1e-5. 
#' @param addD0 a small number that's added for D(X)=0 for numerical evaluation of worst fit (numerical reasons, see details). If addD0=0 the normalized stress for mu!=0 and mu+lambda!=0 is correct, but will give useless normalized stress for mu=0 or mu+lambda!=0.
#' @param principal If 'TRUE', principal axis transformation is applied to the final configuration
#' @param normconf normalize the configuration to sum(delta^2)=1 (as in the power stresses). Default is FALSE. Note that then the distances in confdist do not match manually calculated ones.
#' 
#' @details For numerical reasons with certain parameter combinations, the normalized stress uses a configuration as worst result where every d(X) is 0+addD0. The same number is not added to the delta so there is a small inaccuracy of the normalized stress (but negligible if min(delta)>>addD0). Also, for mu<0 or mu+lambda<0 the normalization cannot generally be trusted (in the worst case of D(X)=0 one would have an 0^(-a)).    
#'
#'
#' @return an object of class 'bcmds' (also inherits from 'smacofP'). It is a list with the components
#' \itemize{
#' \item delta: Observed, untransformed dissimilarities
#' \item tdelta: Observed explicitly transformed dissimilarities, normalized
#' \item dhat: Explicitly transformed dissimilarities (dhats)
#' \item confdist: Configuration dissimilarities
#' \item conf: Matrix of fitted configuration
#' \item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
#' \item ndim: Number of dimensions
#' \item model: Name of MDS model
#' \item type: Must be "ratio" here. 
#' \item niter: Number of iterations
#' \item nobj: Number of objects
#' \item pars: hyperparameter vector theta
#' \item weightmat: 1-diagonal matrix. For compatibility with smacofP classes. 
#' \item parameters, pars, theta: The parameters supplied
#' \item call the call
#' }
#' and some additional components
#' \itemize{
#' \item stress.m: default stress is the explicitly normalized stress on the normalized, transformed dissimilarities
#' \item mu: mu parameter (for attraction)
#' \item lambda: lambda parameter (for repulsion)
#' \item rho: rho parameter (for weights) 
#' }
#'
#' @importFrom stats as.dist dist
#' 
#' @author Lisha Chen & Thomas Rusch
#' 
#' @examples
#' dis<-smacof::kinshipdelta
#' res<-bcmds(dis,mu=2,lambda=1.5,rho=0)
#' res
#' summary(res)
#' plot(res)
#' 
#' @export
bcmds <- function(delta,mu=1,lambda=1,rho=0,type="ratio", ndim=2,weightmat=1-diag(nrow(delta)),itmax=2000,init=NULL,verbose=0,addD0=1e-4,principal=FALSE,normconf=FALSE,acc=1e-5)
{
  if(inherits(delta,"dist") || is.data.frame(delta)) delta <- as.matrix(delta)
  if(!isSymmetric(delta)) stop("Delta is not symmetric.\n")
  if(verbose>0) cat("Minimizing bcStress with mu=",mu,"lambda=",lambda,"rho=",rho,"\n")
  nu <- rho  
  Dorig <- Do <- delta
  n <- nrow(Do)
  if (ndim > (n - 1)) stop("Maximum number of dimensions is n-1!")
  if(is.null(rownames(delta))) rownames(delta) <- 1:n
  labos <- rownames(delta) 
  d <- ndim
  X1 <- init
  xstart <- init
  niter <- itmax
  if(lambda<=0) stop("The lambda parameter must be strictly positive.")
  lambdaorig <- lambda
  lambda <- 1/lambda

  Dnu <- Do^nu
  Dnulam <- Do^(nu+1/lambda)
  diag(Dnu) <- 0
  diag(Dnulam) <- 0
  t <- 0 #not needed
  Grad <- matrix (0, nrow=n, ncol= d)
   if(is.null(X1))
    {
      cmd <- cmds(Do)  
      X1 <- cmd$vec[,1:d]%*%diag(cmd$val[1:d])+enorm(Do)/n/n*0.01*matrix(stats::rnorm(n*d),nrow=n,ncol=d)
      xstart <- X1
    }
  D1 <- as.matrix(dist(X1))
  
  X1 <- X1*enorm(Do)/enorm(D1)

  s1 <- Inf #stressinit
  s0 <- 2 #stress check for update
  stepsize <-0.1
  i <- 0

while ( stepsize > acc && i < niter)
  {
    if (s1 >= s0 && i>1)
    {
       #stepsize if stress (s1) >= old s1 (=s0) oder >2   
       stepsize<- 0.5*stepsize
       X1 <- X0 - stepsize*normgrad
     }
    else 
      {
      #stepsize if stress (s1) < old s1       
      stepsize <- 1.05*stepsize
      X0 <- X1
      D1mu2 <- D1^(mu-2)
      diag(D1mu2) <- 0
      D1mulam2 <- D1^(mu+1/lambda-2)
      diag(D1mulam2) <- 0
      M <- Dnu*D1mulam2-D1mu2*Dnulam    
      E <- matrix(rep(1,n*d),n,d)
      Grad <- X0*(M%*%E)-M%*%X0
      normgrad <- (enorm(X0)/enorm(Grad))*Grad         
      X1 <- X0 - stepsize*normgrad
     }
    i <- i+1
    s0 <- s1 
    D1 <- as.matrix(dist(X1))
    D1mulam <- D1^(mu+1/lambda)
    diag(D1mulam) <- 0
    D1mu <- D1^mu
    diag(D1mu) <- 0
    if(mu+1/lambda==0)
      {
       diag(D1)<-1
       s1 <- sum(Dnu*log(D1))-sum((D1mu-1)*Dnulam)/mu 
      }

    if(mu==0)
      {
      diag(D1)<-1
      s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda) -sum(log(D1)*Dnulam)
      }
    if(mu!=0&(mu+1/lambda)!=0)
    {
        s1 <- sum(Dnu*(D1mulam-1))/(mu+1/lambda)-sum((D1mu-1)*Dnulam)/mu     
    }
    ## Printing and Plotting
    if(verbose > 3 & (i+1)%%100/verbose==0)
      {
        print (paste("niter=",i+1," stress=",round(s1,5), sep=""))
      }

  }
  #For normalization of stress
                                        #Next steps normalize the X1 so that D^*(X) to be smaller than Delta^*
  
  #if((mu!=0) & (mu+lambda)!=0)
  #{
  #    whichd1 <- sum(Dnu*((D1mulam-1)/(mu+1/lambda))^2)
  #    whichd2 <- sum(Dnu*((D1mu-1)/mu)^2)
  #    whichd <- max(whichd1,whichd2)
  #}
  #if(mu==0)
  #{
  #    diag(D1) <- 1
  #    whichd1 <- sum(Dnu*((D1mulam-1)/(mu+lambda))^2)
  #    whichd2 <- sum(Dnu*log(D1)^2)
  #    whichd <- max(whichd1,whichd2)
  #}
  #if(mu+lambda==0)
  #{
  #    diag(D1) <- 1
  #    whichd1 <- sum(Dnu*log(D1)^2)
  #    whichd2 <- sum(Dnu*((D1mu-1)/mu)^2)
  #    whichd <- max(whichd1,whichd2)               
  #}
  #whichd <- sum(Dnu*D1^2)
  #Dlam2 <- Do^(2*1/lambda)
  #diag(Dlam2) <- 0
  #X1 <- X1*sqrt(sum(Dnu*Dlam2))/sqrt(whichd)
  X1 <- X1*sum(Dnu*Do*D1)/sum(Dnu*D1^2)
  D1 <- as.matrix(dist(X1)) 
  D0 <- D1*0+addD0 #for numerical reasons for mu=0 and mu+lambda=0 all get an extra p for "plus"
  Dop <- Do # we could also add addD0 here for reasons of comparability with the D0
  diag(Dop) <- 0
  diag(D0) <- 0
  Dopmulam <- Dop^(mu+1/lambda) #new
  D0mulam <- D0^(mu+1/lambda) #new
  diag(Dopmulam) <- 0
  diag(D0mulam) <- 0
  Dopmu <- Dop^mu #new
  D0mu <- D0^mu #new
  diag(Dopmu) <- 0 #new
  diag(D0mu) <- 0 #new
  Dpnu <-  Dop^nu
  Dpnulam <- Dop^(nu+1/lambda)
  diag(Dpnu) <- 0
  diag(Dpnulam) <- 0
  if(mu+1/lambda==0) {
       diag(D0) <- 1
       diag(Do)<-1 #new
       norm0 <- sum(Dnu*log(D0))-sum((D0mu-1)*Dnulam)/mu
       normo <- sum(Dnu*log(Do))-sum((Dopmu-1)*Dnulam)/mu
       s1n <- (s1-normo)/(norm0-normo)       
       }
  if(mu==0) {
       diag(D0) <- 1
       diag(Do)<-1 #new
       norm0 <- sum(Dnu*(D0mulam-1))/(mu+1/lambda) - sum(log(D0)*Dnulam)
       normo <- sum(Dnu*(Dopmulam-1))/(mu+1/lambda) - sum(log(Do)*Dnulam)
      s1n <- (s1-normo)/(norm0-normo)
      }
  if(mu!=0&(mu+1/lambda)!=0)
  {
      #D0 <- D1*0
      #D0mu <- D0^mu # infinity
      #diag(D0mu) <- 0 
      #D0mulam <- D0^(mu+1/lambda) 
      #diag(D0mulam) <- 0
      #Domulam <- Do^(mu+1/lambda) 
      #diag(Domulam) <- 0
      #Domu <- Do^mu #new
      #diag(Domu) <- 0 #new
      normo <- sum(Dpnu*(Dopmulam-1))/(mu+1/lambda)-sum((Dopmu-1)*Dpnulam)/mu
      norm0 <- sum(Dpnu*(D0mulam-1))/(mu+1/lambda)-sum((D0mu-1)*Dpnulam)/mu      
      #norm0 <- sum(Dnu*(D0mulam-1))/(mu+1/lambda)-sum((D0mu-1)*Dnulam)/mu
      #normo <- sum(Dnu*(Domulam-1))/(mu+1/lambda)-sum((Domu-1)*Dnulam)/mu
      s1n <- (s1-normo)/(norm0-normo)
      # s1n <- 1-s1/normo #normalized stress
  }
  result <- list()
  result$delta <- stats::as.dist(Dorig)
  result$dhat <- stats::as.dist(Do)  #TODO: Check again
  if(isTRUE(normconf)) X1 <- X1/enorm(X1)
  if (principal) {
        X1_svd <- svd(X1)
        X1 <- X1 %*% X1_svd$v
  }
  attr(X1,"dimnames")[[1]] <- labos
  attr(X1,"dimnames")[[2]] <- paste("D",1:ndim,sep="")  
  result$iord <- order(as.vector(Do)) ##TODO: Check again
  result$confdist <- stats::as.dist(D1) #TODO: Check again
  result$conf <- X1 #new
  result$stress <- sqrt(s1n)
  weightmat <- stats::as.dist(1-diag(n))
  spoint <- spp(result$delta, result$confdist, weightmat)
  resmat<-spoint$resmat
  rss <- sum(spoint$resmat[lower.tri(spoint$resmat)])
  spp <- spoint$spp
  result$spp <- spp
  result$ndim <- ndim
  result$weightmat <- weightmat
  result$resmat <- resmat
  result$rss <- rss
  result$init <- xstart
  result$model<- "Box-Cox MDS"
  result$niter <- i
  result$nobj <- n
  result$type <- "ratio"
  result$call <- match.call()
  result$stress.m <- s1n
  result$stress.r <- s1
  result$tdelta <- stats::as.dist(Do)
  result$parameters <- c(mu=mu,lambda=lambdaorig,rho=nu)
  result$pars <- c(mu=mu,lambda=lambdaorig,rho=nu)
  result$theta <- c(mu=mu,lambda=lambdaorig,rho=nu)
  #result$tweightmat <- weightmat
  result$mu <- mu
  result$lambda <- lambdaorig
  result$rho <- nu
  class(result) <- c("bcmds","smacofP","smacofB","smacof")
  return(result)
}


#' Classical Scaling
#' @param Do dissimilarity matrix
cmds <- function(Do)
  {
    n <- nrow(Do)
    J <- diag(rep(1,n)) - 1/n * rep(1,n) %*% t(rep(1,n))
    B <- - 1/2 * J %*% (Do^2) %*% J
    pc <- eigen(B)
    return(pc)
  }


#' @rdname bcmds
bcStressMin <- bcmds

#' @rdname bcmds
bcstressMin <- bcmds

#' @rdname bcmds
boxcoxmds <- bcmds



## function (delta, init = NULL, verbose = 0, ndim = 2, lambda = 1, 
##     mu = 1, nu = 0, itmax = 10000) 
## {
##     Do <- delta
##     d <- ndim
##     X1 <- init
##     niter <- itmax
##     lambdaorig <- lambda
##     lambda <- 1/lambda
##     n <- nrow(Do)
##     Dnu <- Do^nu
##     Dnulam <- Do^(nu + 1/lambda)
##     diag(Dnu) <- 0
##     diag(Dnulam) <- 0
##     t <- 0
##     Grad <- matrix(0, nrow = n, ncol = d)
##     if (is.null(X1)) {
##         cmd <- cmds(Do)
##         X1 <- cmd$vec[, 1:d] %*% diag(cmd$val[1:d]) + enorm(Do)/n/n * 
##             0.01 * matrix(stats::rnorm(n * d), nrow = n, ncol = d)
##     }
##     D1 <- as.matrix(dist(X1))
##     X1 <- X1 * enorm(Do)/enorm(D1)
##     s1 <- Inf
##     s0 <- 2
##     stepsize <- 0.1
##     i <- 0
##     while (stepsize > 1e-05 && i < niter) {
##         if (s1 >= s0 && i > 1) {
##             stepsize <- 0.5 * stepsize
##             X1 <- X0 - stepsize * normgrad
##         }
##         else {
##             stepsize <- 1.05 * stepsize
##             X0 <- X1
##             D1mu2 <- D1^(mu - 2)
##             diag(D1mu2) <- 0
##             D1mulam2 <- D1^(mu + 1/lambda - 2)
##             diag(D1mulam2) <- 0
##             M <- Dnu * D1mulam2 - D1mu2 * Dnulam
##             E <- matrix(rep(1, n * d), n, d)
##             Grad <- X0 * (M %*% E) - M %*% X0
##             normgrad <- (enorm(X0)/enorm(Grad)) * Grad
##             X1 <- X0 - stepsize * normgrad
##         }
##         i <- i + 1
##         s0 <- s1
##         D1 <- as.matrix(dist(X1))
##         D1mulam <- D1^(mu + 1/lambda)
##         diag(D1mulam) <- 0
##         D1mu <- D1^mu
##         diag(D1mu) <- 0
##         if (mu + 1/lambda == 0) {
##             diag(D1) <- 1
##             s1 <- sum(Dnu * log(D1)) - sum((D1mu - 1) * Dnulam)/mu
##         }
##         if (mu == 0) {
##             diag(D1) <- 1
##             s1 <- sum(Dnu * (D1mulam - 1))/(mu + 1/lambda) - 
##                 sum(log(D1) * Dnulam)
##         }
##         if (mu != 0 & (mu + 1/lambda) != 0) {
##             s1 <- sum(Dnu * (D1mulam - 1))/(mu + 1/lambda) - 
##                 sum((D1mu - 1) * Dnulam)/mu
##         }
##         if (verbose > 3 & (i + 1)%%100/verbose == 0) {
##             print(paste("niter=", i + 1, " stress=", round(s1, 
##                 5), sep = ""))
##         }
##     }
##    # X1 <- X1 * sum(Do * D1)/sum(D1^2)
##     D1 <- as.matrix(dist(X1))
##     D0 <- D1 * 0 + 1e-04
##     Dop <- Do + 1e-04
##     Dopmulam <- Dop^(mu + 1/lambda)
##     Dopmu <- Dop^mu
##     diag(Dopmu) <- 0
##     diag(Dop) <- 1
##     diag(D0) <- 1
##     D0mu <- D0^mu
##     Dopmu <- Dop^mu
##     diag(Dopmu) <- 0
##     Dpnu <- Dop^nu
##     Dpnulam <- Dop^(nu + 1/lambda)
##     if (mu + 1/lambda == 0) {
##         norm0 <- sum(Dpnu * log(D0)) - sum((D0mu - 1) * Dpnulam)/mu
##         normo <- sum(Dpnu * log(Dop)) - sum((Dopmu - 1) * Dpnulam)/mu
##         s1n <- (s1 - normo)/(norm0 - normo)
##     }
##     if (mu == 0) {
##         norm0 <- sum(Dpnu * (D0mulam - 1))/(mu + 1/lambda) - 
##             sum(log(D0) * Dpnulam)
##         normo <- sum(Dpnu * (Dopmulam - 1))/(mu + 1/lambda) - 
##             sum(log(Dop) * Dpnulam)
##         s1n <- (s1 - normo)/(norm0 - normo)
##     }
##     if (mu != 0 & (mu + 1/lambda) != 0) {
##         norm0 <- sum(Dpnu * (D0mulam - 1))/(mu + 1/lambda) - 
##             sum((D0mu - 1) * Dpnulam)/mu
##         normo <- sum(Dpnu * (Dopmulam - 1))/(mu + 1/lambda) - 
##             sum((Dopmu - 1) * Dpnulam)/mu
##         s1n <- (s1 - normo)/(norm0 - normo)
##     }
##     result <- list()
##     result$conf <- X1
##     result$confdist <- stats::as.dist(D1)
##     result$delta <- stats::as.dist(Do)
##     result$obsdiss <- stats::as.dist(Do)
##     result$mu <- mu
##     result$lambda <- lambda
##     result$nu <- nu
##     result$pars <- c(mu, lambda, nu)
##     result$model <- "Box-Cox Stress MDS"
##     result$call <- match.call()
##     result$ndim <- ndim
##     result$nobj <- n
##     result$niter <- i
##     result$theta <- c(mu, lambda, nu)
##     result$stress.r <- s1
##     result$stress.m <- s1n
##     result$stress <- sqrt(s1n)
##     result$type <- "Box-Cox Stress"
##     class(result) <- c("smacofP", "smacofB", "smacof")
##     return(result)
## }
