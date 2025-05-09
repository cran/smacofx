% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmds.R
\name{lmds}
\alias{lmds}
\title{Local MDS}
\usage{
lmds(
  delta,
  k = 2,
  tau = 1,
  type = "ratio",
  ndim = 2,
  weightmat = 1 - diag(nrow(delta)),
  itmax = 5000,
  acc = 1e-05,
  init = NULL,
  verbose = 0,
  principal = FALSE,
  normconf = FALSE
)
}
\arguments{
\item{delta}{dissimilarity or distance matrix, dissimilarity or distance data frame or 'dist' object}

\item{k}{the k neighbourhood parameter}

\item{tau}{the penalty parameter (suggested to be in [0,1])}

\item{type}{what type of MDS to fit. Only "ratio" currently.}

\item{ndim}{the dimension of the configuration}

\item{weightmat}{a matrix of finite weights. Not implemented.}

\item{itmax}{number of optimizing iterations, defaults to 5000.}

\item{acc}{accuracy (lowest stepsize). Defaults to 1e-5.}

\item{init}{initial configuration. If NULL a classical scaling solution is used.}

\item{verbose}{prints info if > 0 and progress if > 1.}

\item{principal}{If 'TRUE', principal axis transformation is applied to the final configuration}

\item{normconf}{normalize the configuration to sum(delta^2)=1 (as in the power stresses). Note that then the distances in confdist do not match the manually calculated ones.}
}
\value{
an object of class 'lmds' (also inherits from 'smacofP'). See \code{\link{powerStressMin}}. It is a list with the components as in power stress
\itemize{
\item delta: Observed, untransformed dissimilarities
\item tdelta: Observed explicitly transformed dissimilarities, normalized
\item dhat: Explicitly transformed dissimilarities (dhats)
\item confdist: Configuration dissimilarities
\item conf: Matrix of fitted configuration
\item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
\item ndim: Number of dimensions
\item model: Name of MDS model
\item type: Is "ratio" here. 
\item niter: Number of iterations
\item nobj: Number of objects
\item pars: explicit transformations hyperparameter vector theta
\item weightmat: 1-diagonal matrix (for compatibility with smacof classes)
\item parameters, pars, theta: The parameters supplied
\item call the call
}
and some additional components
\itemize{
\item stress.m: default stress is the explicitly normalized stress on the normalized, transformed dissimilarities
\item tau: tau parameter
\item k: k parameter
}
}
\description{
This function minimizes the Local MDS Stress of Chen & Buja (2006) via gradient descent. This is a ratio metric scaling method.
}
\details{
Note that k and tau are not independent. It is possible for normalized stress to become negative if the tau and k combination is so that the absolute repulsion for the found configuration dominates the local stress substantially less than the repulsion term does for the solution of D(X)=Delta, so that the local stress difference between the found solution and perfect solution is nullified. This can typically be avoided if tau is between 0 and 1. If not, set k and or tau to a smaller value.
}
\examples{
dis<-smacof::kinshipdelta
res<- lmds(dis,k=2,tau=0.1)
res
summary(res)
plot(res)

}
\author{
Lisha Chen & Thomas Rusch
}
