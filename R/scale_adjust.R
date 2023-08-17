#' Adjusts a configuration
#' 
#'@param conf a configuration
#'@param ref a reference configuration (only for scale="proc") 
#'@param scale Scale adjustment. "std" standardizes each column of the configurations to mean=0 and sd=1, "sd" scales the configuration by the maximum standard devation of any column, "proc" adjusts the fitted configuration to the reference
#' @return The scale adjusted configuration.
#'
#' @importFrom stats sd
#' @importFrom smacof Procrustes 
scale_adjust <- function(conf,ref,scale=c("sd","std","proc","none"))
{
  #conf target, ref reference
  if(scale=="std") conf <- base::scale(conf)
  if(scale=="sd") conf <- conf/max(apply(conf,2,stats::sd))
  if(scale=="proc") conf <- smacof::Procrustes(ref,conf)$Yhat
  if(scale=="none") conf <- conf
  conf
}
