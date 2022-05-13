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
#' potential=function(x){sum(0.5*x^2/sigma^2)}
#' grad=function(x){x/sigma^2}
#' n=10^4
#' g=1.5
#' h=0.20
#' T_phys=2
#' x_init=rnorm(d)*sigma
#' output=malt(potential,grad,n,g,h,T_phys,x_init)
#' ess(output)
ess=function(output){
  ess_mean=coda::effectiveSize(output$chain)
  ess_square=coda::effectiveSize(output$chain^2)
  return(list(ess_mean=ess_mean,ess_square=ess_square))
}
