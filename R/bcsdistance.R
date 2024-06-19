#' Calculates the blended Chi-square distance matrix between n vectors
#'
#' The pairwise blended chi-distance of two vectors x and y is sqrt(sum(((x[i]-y[i])^2)/(2*(ax[i]+by[i])))), with originally a in [0,1] and b=1-a as in Lindsay (1994) (but we allow any non-negative a and b). The function calculates this for all pairs of rows of a matrix or data frame x.    
#'
#' @param x an n times p numeric matrix or data frame. Note that the valeus of x must be non-negative.
#' @param a first blending weight. Must be non-negative and should be in [0,1] if a blended chi-square distance as in Lindsay (1994) is sought. Defaults to 0.5.
#' @param b second blending weight. Must be non-negative and should be 1-a if a blended chi-square distance as in Lindsay (1994) is sought. Defaults to 1-a.
#'
#' @return a symmetric n times n matrix of pairwise blended chi-square distance (between rows of x) with 0 in the main diagonal. It is an object of class distance and matrix with attributes "method", "type" and "par", the latter returning the a and b values.  
#'
#' @references Lindsay (1994). Efficiency versus robustness: the case for minimum Hellinger distance and related methods. Annals of Statistics, 22 (2), 1081-1114. <doi:10.1214/aos/1176325512>
#' 
#' @export
bcsdistance <- function(x,a=0.5,b=1-a)
{
  #inspired by oldDistance in package analogue.
  if(any(x <0)) stop("Blended Chi-Square Distance can only be calculated for non-negative values in x.")
  if(a < 0 || b < 0) stop("Blended Chi-Square Distance can only be calculated for non-negative values of a and b.")
  b.chisq <- function(x, y, a, b)
      {
        inds <- !(x == 0L & y == 0L)
        sqrt(sum(((x[inds] - y[inds])^2) / (2*(a*x[inds] + b*y[inds]))))
      }
  b.chi.square <- function(y,x,a,b) apply(x, 1, b.chisq, y, a, b)
  y <- x
  n.vars <- ncol(x)
  facs.x <- facs.y <- rep(FALSE, n.vars)
  x.names <- rownames(x)
  x <- data.matrix(x)
  y.names <- rownames(y)
  y <- data.matrix(y)
  dimx <- dim(x)
  dimy <- dim(y)
  dimnames(x) <- dimnames(y) <- NULL
  res <-  apply(y, 1, b.chi.square, x, a, b)
  if(is.null(dim(res))) {
      names(res) <- x.names
   } else {
      colnames(res) <- y.names
      rownames(res) <- x.names
    }
  attr(res, "method") <- "blended chi-square"
  attr(res, "type") <- "symmetric"
  attr(res, "par") <- c("a"=a,"b"=b)
  class(res) <- c("distance","matrix")
  return(res)
 }
