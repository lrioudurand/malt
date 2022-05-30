#' ess
#'
#' @description computes effective sample sizes of the chain
#'
#' @param output an output of the malt function
#'
#' @return a list of effective sample sizes when estimating the first and second moments
#' @importFrom coda effectiveSize
#' @export
#'
#' @examples
#' d=50
#' sigma=((d:1)/d)^(1/2)
#' init=rnorm(d)*sigma
#' U=function(x){sum(0.5*x^2/sigma^2)}
#' grad=function(x){x/sigma^2}
#' n=10^4
#' g=1.5
#' h=0.20
#' L=10
#' output=malt(init,U,grad,n,g,h,L)
#' ess(output)
ess=function(output){
  ess_mean=coda::effectiveSize(output$samples)
  ess_square=coda::effectiveSize(output$samples^2)
  return(list(ess_mean=ess_mean,ess_square=ess_square))
}
