% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.smacofP.R
\name{plot.smacofP}
\alias{plot.smacofP}
\title{S3 plot method for smacofP objects}
\usage{
\method{plot}{smacofP}(
  x,
  plot.type = "confplot",
  plot.dim = c(1, 2),
  bubscale = 1,
  col,
  label.conf = list(label = TRUE, pos = 3, col = 1, cex = 0.8),
  hull.conf = list(hull = FALSE, col = 1, lwd = 1, ind = NULL),
  shepard.x = NULL,
  identify = FALSE,
  type = "p",
  cex = 0.5,
  pch = 20,
  asp = 1,
  main,
  xlab,
  ylab,
  xlim,
  ylim,
  col.hist = NULL,
  legend = TRUE,
  legpos,
  loess = TRUE,
  shepard.lin = TRUE,
  ...
)
}
\arguments{
\item{x}{an object of class smacofP}

\item{plot.type}{String indicating which type of plot to be produced: "confplot", "resplot", "Shepard", "stressplot","transplot", "bubbleplot" (see details)}

\item{plot.dim}{dimensions to be plotted in confplot; defaults to c(1, 2)}

\item{bubscale}{Scaling factor (size) for the bubble plot}

\item{col}{vector of colors for the points}

\item{label.conf}{List with arguments for plotting the labels of the configurations in a configuration plot (logical value whether to plot labels or not, label position, label color)}

\item{hull.conf}{Option to add convex hulls to a configuration plot. Hull index needs to be provided.}

\item{shepard.x}{Shepard plot only: original data (e.g. correlation matrix) can be provided for plotting on x-axis}

\item{identify}{If 'TRUE', the 'identify()' function is called internally that allows to add configuration labels by mouse click}

\item{type}{What type of plot should be drawn (see also 'plot')}

\item{cex}{Symbol size.}

\item{pch}{Plot symbol}

\item{asp}{Aspect ratio; defaults to 1 so distances between x and y are represented accurately; can lead to slighlty weird looking plots if the variance on one axis is much smaller than on the other axis; use NA if the standard type of R plot is wanted where the ylim and xlim arguments define the aspect ratio - but then the distances seen are no longer accurate}

\item{main}{plot title}

\item{xlab}{label of x axis}

\item{ylab}{label of y axis}

\item{xlim}{scale of x axis}

\item{ylim}{scale of y axis}

\item{col.hist}{Color of the borders of the histogram.}

\item{legend}{Flag whether legends should be drawn for plots that have legends}

\item{legpos}{Position of legend in plots with legends}

\item{loess}{if TRUE a loess fit (by Tukey's rescending M-Estimator) of configuration distances explained by delta  is added to the Shepard plot}

\item{shepard.lin}{Shepard plot only: if TRUE the Shepard plot is linearized so d^kappa~delta^lambda. If FALSE d~delta^lambda}

\item{...}{Further plot arguments passed: see 'plot.smacof' and 'plot' for detailed information.}
}
\value{
no return value; just plots for class 'smacofP' (see details)
}
\description{
S3 plot method for smacofP objects
}
\details{
\itemize{
\item  Configuration plot (plot.type = "confplot"): Plots the MDS configuration.
 \item Residual plot (plot.type = "resplot"): Plots the dhats f(T(delta)) against the transformed fitted distances T(d(X)).
 \item (Linearized) Shepard diagram (plot.type = "Shepard"): Is shep.lin=TRUE a diagram with the transformed observed normalized dissimilarities (T(delta) on x)  against the transformed fitted distance (T(d(X) on y) as well as a loess curve and a regression line corresponding to type (linear without intercept for ratio, linear for interval and isotonic for ordinal). If shep.lin=FALSE it uses the untransformed delta. Note that the regression line corresponds to the optimal scaling results (dhat) only up to a linear transformation. 
 \item Transformation Plot (plot.type = "transplot"): Diagram with normalized observed dissimilarities (delta, light grey) and the normalized explicitly transformed dissimilarities (T(Delta), darker) against the untransformed fitted distances (d(X)) together with a nonlinear regression curve corresponding to the explicit transformation (fitted power transformation). This is most useful for ratio models with power transformations as the transformations can be read of directly. For other MDS models and stresses, it still gives a quick way to assess how the explicit transformations worked.  
 \item Stress decomposition plot (plot.type = "stressplot"): Plots the stress contribution in of each observation. Note that it rescales the stress-per-point (SPP) from the corresponding function to percentages (sum is 100). The higher the contribution, the worse the fit.
 \item Bubble plot (plot.type = "bubbleplot"): Combines the configuration plot with the point stress contribution. The larger the bubbles, the worse the fit.
\item histogram (‘plot.type = "histogram"’: gives a weighted histogram of the dissimilarities (weighted with tweightmat if exists else with weightmat). For optional arguments, see ‘wtd.hist’.
}
}
\examples{
dis<-as.matrix(smacof::kinshipdelta)
res<-powerStressMin(dis)
plot(res)
plot(res,"Shepard")
plot(res,"resplot")
plot(res,"transplot")
plot(res,"stressplot")
plot(res,"bubbleplot")
plot(res,"histogram")
}
