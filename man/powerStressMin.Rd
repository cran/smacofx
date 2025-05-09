% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/powerStressMin.R
\name{powerStressMin}
\alias{powerStressMin}
\alias{powerstressMin}
\alias{postmds}
\alias{pstressMin}
\alias{pStressMin}
\alias{pstressmds}
\title{Power Stress SMACOF}
\usage{
powerStressMin(
  delta,
  kappa = 1,
  lambda = 1,
  nu = 1,
  type = "ratio",
  weightmat = 1 - diag(nrow(delta)),
  init = NULL,
  ndim = 2,
  acc = 1e-06,
  itmax = 10000,
  verbose = FALSE,
  principal = FALSE
)

powerstressMin(
  delta,
  kappa = 1,
  lambda = 1,
  nu = 1,
  type = "ratio",
  weightmat = 1 - diag(nrow(delta)),
  init = NULL,
  ndim = 2,
  acc = 1e-06,
  itmax = 10000,
  verbose = FALSE,
  principal = FALSE
)

postmds(
  delta,
  kappa = 1,
  lambda = 1,
  nu = 1,
  type = "ratio",
  weightmat = 1 - diag(nrow(delta)),
  init = NULL,
  ndim = 2,
  acc = 1e-06,
  itmax = 10000,
  verbose = FALSE,
  principal = FALSE
)

pstressMin(
  delta,
  kappa = 1,
  lambda = 1,
  nu = 1,
  type = "ratio",
  weightmat = 1 - diag(nrow(delta)),
  init = NULL,
  ndim = 2,
  acc = 1e-06,
  itmax = 10000,
  verbose = FALSE,
  principal = FALSE
)

pStressMin(
  delta,
  kappa = 1,
  lambda = 1,
  nu = 1,
  type = "ratio",
  weightmat = 1 - diag(nrow(delta)),
  init = NULL,
  ndim = 2,
  acc = 1e-06,
  itmax = 10000,
  verbose = FALSE,
  principal = FALSE
)

pstressmds(
  delta,
  kappa = 1,
  lambda = 1,
  nu = 1,
  type = "ratio",
  weightmat = 1 - diag(nrow(delta)),
  init = NULL,
  ndim = 2,
  acc = 1e-06,
  itmax = 10000,
  verbose = FALSE,
  principal = FALSE
)
}
\arguments{
\item{delta}{dist object or a symmetric, numeric data.frame or matrix of distances}

\item{kappa}{power of the transformation of the fitted distances; defaults to 1}

\item{lambda}{the power of the transformation of the proximities; defaults to 1}

\item{nu}{the power of the transformation for weightmat; defaults to 1}

\item{type}{what type of MDS to fit. One of "ratio" or "interval". Default is "ratio".}

\item{weightmat}{a matrix of finite weights or dist object}

\item{init}{starting configuration}

\item{ndim}{dimension of the configuration; defaults to 2}

\item{acc}{numeric accuracy of the iteration. Default is 1e-6.}

\item{itmax}{maximum number of iterations. Default is 10000.}

\item{verbose}{should internal messages be printed; if > 0 then yes (iteration progress with >2)}

\item{principal}{If 'TRUE', principal axis transformation is applied to the final configuration}
}
\value{
a 'smacofP' object (inheriting from 'smacofB', see \code{\link[smacof]{smacofSym}}). It is a list with the components
\itemize{
\item delta: Observed, untransformed dissimilarities
\item tdelta: Observed explicitly transformed dissimilarities, normalized
\item dhat: Explicitly transformed dissimilarities (dhats), optimally scaled and normalized 
\item confdist: Transformed fitted configuration distances
\item conf: Matrix of fitted configuration
\item stress: Default stress  (stress 1; sqrt of explicitly normalized stress)
\item spp: Stress per point 
\item ndim: Number of dimensions
\item model: Name of smacof model
\item niter: Number of iterations
\item nobj: Number of objects
\item type: Type of MDS model 
\item weightmat: weighting matrix as supplied 
\item stress.m: Default stress (stress-1^2)
\item tweightmat: transformed weighthingmatrix (here weightmat^nu)
}
}
\description{
An implementation to minimize power stress by majorization with ratio or interval optimal scaling. Usually more accurate but slower than powerStressFast. Uses a repeat loop.
}
\examples{
dis<-smacof::kinshipdelta
res<-powerStressMin(dis,type="ratio",kappa=2,lambda=1.5,itmax=1000)
res
summary(res)
plot(res)

}
\seealso{
\code{\link[smacof]{smacofSym}}
}
