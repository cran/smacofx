#' smacofx: Flexible multidimensional scaling methods and SMACOF extensions 
#' 
#' Flexible multidimensional scaling (MDS) methods centered around the Majorization algorithm. The package contains various functions, wrappers, methods and classes for fitting, plotting and displaying a large number of different flexible MDS models such as Torgerson scaling, ratio, interval and nonmetric MDS with majorization, nonlinear MDS with optimal powers for dissimilarities, Sammon mapping with ratio and interval optimal scaling, multiscale MDS with ratio and interval optimal scaling, Alscal (s-stress) MDS with ratio and interval optimal scaling, elastic scaling with ratio and interval optimal scaling, r-stress MDS for ratio, interval and nonmetric scaling, power stress for interval and ratio optimal scaling, restricted power-stress with ratio and interval optimal scaling, approximate power-stress with ratio scaling, curvilinear component analysis with ratio, interval and ordinal optimal scaling, power curvilinear component analysis with ratio, interval and ordinal optimal scaling, Box-Cox MDS and local MDS. Some functions are suitably flexible to allow any other sensible combination of explicit power transformations for weights, distances and input proximities with implicit ratio, interval or ordinal optimal scaling of the input proximities. Most functions use a majorization algorithm.
#'
#' The package provides:
#'
#' Models:
#' \itemize{
#' \item alscal... ALSCAL (s-stress) MDS with ratio, interval optimal scaling
#' \item elscal.. Elastic scaling MDS with ratio, interval optimal scaling
#' \item multiscale... Multiscale MDSwith ratio, interval optimal scaling
#' \item rstressMin .. R-Stress MDS with ratio, interval, ordinal optimal scaling
#' \item powerStressMin... power stress MDS (POST-MDS) with ratio, interval optimal scaling
#' \item apStressMin... approximate POST-MDS with ratio, interval optimal scaling
#' \item rpowerStressMin... restricted POST-MDS with ratio, interval optimal scaling
#' \item clca ... curvilinear component analysis with ratio optimal scaling
#' \item clda ... curvilinear distance analysis with ratio optimal scaling (ie.e. clca with Isomap distances)
#' \item bcmds ... Box-Cox MDS with ratio optimal scaling
#' \item lmds... Local MDS with ratio optimal scaling
#' \item sammonmap... Sammon mapping with ratio, interval optimal scaling
#' \item smds ... power curvilinear component analysis with ratio, interval, ordinal optimal scaling
#' \item pclca ... power curvilinear component analysis with ratio, interval, ordinal optimal scaling
#' \item smdda ... sparsified multidimensional distance analsysis with ratio interval scaling (this is smds with Isomap distances) 
#' \item spmdda ... sparsified power multidimensional distance analsysis with ratio interval scaling (this is spmds with Isomap distances)  
#' \item smds ... sparsified multidimensional scaling with ratio interval scaling (inspired by CLCA, so fitted distances larger than tau are weighted with 0)
#' \item spmds ... sparsified multidimensional scaling with ratio interval scaling (inspired by CLCA, so fitted distances larger than tau are weighted with 0)
#' \item opmds ... nonlinear MDS with optimal power for the dissimilarities.
#' }
#'
#' 
#' Classes and Methods: 
#' The objects are of classes that extend the S3 classes smacof and smacofB. For the objects returned by the high-level functions S3 methods for standard generics were implemented, including print, coef, residuals, summary, plot, plot3dstatic.

#' Wrappers and convenience functions for the model objects:
#' \itemize{
#' \item bootmds ... Bootstrapping an MDS model
#' \item biplotmds ... MDS Biplots
#' \item icExploreGen ... Expore initial configurations 
#' \item jackmds ... Jackknife for MDS 
#' \item multistart ... Multistart function for MDS  
#' \item permtest ... Permutation test for MDS
#'}
#'
#' Wrappers:
#' \itemize{
#' \item cmdscale ... stats::cmdscale but returns an S3 objects to be used with smacof classes 
#' \item sammon... MASS::sammon but returns S3 objects to be used with smacof classes 
#'}
#'
#'
#' @examples
#' \donttest{
#' data(BankingCrisesDistances)
#'
#' res<-rStressMin(BankingCrisesDistances[,1:69],type="ordinal",r=2)
#' res
#' 
#' summary(res)
#' plot(res)
#' plot(res,"transplot")
#' plot(res,"Shepard")
#'
#' msres<- multistart(res)
#'
#' res2<-msres$best
#' permtest(res2)
#'}
#' @keywords internal
#' @aliases smacofx-package 
"_PACKAGE"
NULL
