#' malt
#'
#' @description Generates approximate samples from a distribution with density \deqn{\Pi(x)\propto e^{-U(x)}}{Pi(x)=exp(-U(x))/C} Given a potential function \eqn{U}{U} and its gradient evaluation, implements the sampling algorithm: Metropolis Adjusted Langevin Trajectories (Riou-Durand and Vogrinc 2022). Details available at: https://arxiv.org/abs/2202.13230.
#'
#' @param potential a function (potential evaluation)
#' @param grad a function (gradient evaluation)
#' @param n a positive integer (number of MCMC samples)
#' @param g a non-negative number (friction)
#' @param h a positive number (time step)
#' @param T_phys a positive number (integration time)
#' @param x_init a vector (initial position)
#' @param warm_up a boolean (should the chain be warmed up? default=FALSE). If TRUE then unadjusted trajectories are drawn for n/2 iterations followed by adjusted trajectories for n/2 iterations. The chain starts from the resulting output and the warm up phase is discarded.
#'
#' @return a list with the samples from the chain, the acceptance rate, the effective sample sizes for estimating the first and second moments of each component.
#' @importFrom stats rnorm rexp var
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
#' malt(potential,grad,n,g,h,T_phys,x_init)
malt=function(potential,grad,n,g,h,T_phys,x_init,warm_up=F){
  d=length(x_init)
  if(T_phys<h){warning("Integration time is lower than the time step, setting T_phys=h instead.");T_phys=h}
  L=floor(T_phys/h)
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  chain=matrix(0,nrow=n,ncol=d)
  malt_v=rnorm(d)
  malt_x=x_init
  chain[1,]=malt_x
  #warm up (doubles the computational time): draw unadjusted trajectories for n/2 iterations then adjusted trajectories for n/2 iterations.
  if(warm_up==T){
    for(i in 1:(n-1)){
      x=malt_x
      grad_x=grad(malt_x)
      if(i>(n/2)){
      potential_malt_x=potential(malt_x)
      sum_grad_sq=sum(grad_x^2)
      }
      v=rnorm(d)
      Delta=0
      for(j in 1:L){
        v=v-(h/2)*grad_x
        x=x+h*v
        grad_y=grad(x)
        Delta=Delta-(h/2)*sum(v*(grad_x+grad_y))
        v=v-(h/2)*grad_y
        v=eta*v+zeta*rnorm(d)
        grad_x=grad_y
      }
      if(i<=n/2){malt_x=x}else{
        log_alpha=(h^2/8)*(sum(grad_x^2)-sum_grad_sq)+potential(x)-potential_malt_x+Delta
        if(rexp(1)>log_alpha){malt_x=x}
        }
    }
    chain[1,]=malt_x
  }
  for(i in 1:(n-1)){
    x=malt_x
    grad_x=grad(malt_x)
    potential_malt_x=potential(malt_x)
    sum_grad_sq=sum(grad_x^2)
    v=rnorm(d)
    Delta=0
    for(j in 1:L){
      v=v-(h/2)*grad_x
      x=x+h*v
      grad_y=grad(x)
      Delta=Delta-(h/2)*sum(v*(grad_x+grad_y))
      v=v-(h/2)*grad_y
      v=eta*v+zeta*stats::rnorm(d)
      grad_x=grad_y
    }
    log_alpha=(h^2/8)*(sum(grad_x^2)-sum_grad_sq)+potential(x)-potential_malt_x+Delta
    if(stats::rexp(1)>log_alpha){malt_x=x}
    chain[i+1,]=malt_x
  }
  return(list(chain=chain, acceptance_rate=mean(chain[2:n,1]!=chain[1:(n-1),1]),grad_steps=L,means=apply(chain,2,mean),variances=apply(chain,2,var)))
}
