#' malt
#'
#' @description Generates approximate samples from a distribution with density \deqn{\Pi(x)\propto e^{-U(x)}}{Pi(x)=exp(-U(x))/C} Given a potential function \eqn{U}{U} and its gradient evaluation, implements the sampling algorithm: Metropolis Adjusted Langevin Trajectories (Riou-Durand and Vogrinc 2022). Details available at: https://arxiv.org/abs/2202.13230.
#'
#' @param init Real vector. Initial values for the sampling algorithm.
#' @param U A potential function to return the log-density of the distribution to be sampled from, up to an additive constant. It should input a real vector of the same length as init and output a scalar.
#' @param grad A function to return the gradient of the potential. It should input and output a real vector of the same length as init.
#' @param n Positive integer. The number of samples to be generated.
#' @param g Non-negative real number. The friction, a.k.a damping parameter. The choice g=0 boils down to Hamiltonian Monte Carlo.
#' @param h Positive real number. The time step.
#' @param L Positive integer. The number of steps per trajectory. The choice L=1 boils down to the Metropolis Adjusted Langevin Algorithm.
#' @param warm_up Logical. Should the chain be warmed up? If TRUE, the samples are generated after a warm-up phase of n successive trajectories. The first half of the warm-up phase is composed by unadjusted trajectories.
#'
#' @return Returns a list with the following objects:
#'   \item {samples}{a matrix whose rows are the samples generated.}
#'   \item {draw}{a vector corresponding to the last draw of the chain.}
#'   \item {acceptance}{the acceptance rate of the chain. An acceptance rate close to zero/one indicates that the time step chosen is respectively too large/small. Optimally, the time step should be tuned to obtain an acceptance rate slightly above 65%.}
#' @importFrom stats rnorm rexp
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
#' malt(init,U,grad,n,g,h,L)
malt=function(init,U,grad,n,g,h,L,warm_up=F){
  d=length(init)
  if(d>=2){}else{stop("length(init) should be greater or equal than 2.")}
  if(n>=1){n=floor(n)}else{warning("The number of samples should be a positive integer, setting n=10000 instead.");n=10^4}
  if(g>=0){g=as.numeric(g)}else{warning("The friction/damping should be a non-negative real number, setting g=0 instead.");g=0}
  if(h>0){h=as.numeric(h)}else{stop("The time step should be a positive real number.")}
  if(L>=1){L=floor(L)}else{warning("The number of steps should be a positive integer, setting L=1 instead.");L=1}
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  half=h/2
  small=h^2/8
  chain=matrix(0,nrow=n,ncol=d)
  malt_x=init
  chain[1,]=malt_x
  #warm up (doubles the computational time): draw unadjusted trajectories for n/2 iterations then adjusted trajectories for n/2 iterations.
  if(warm_up==T){
    for(i in 1:(n-1)){
      x=malt_x
      grad_x=grad(malt_x)
      if(i>(n/2)){
      U_malt_x=U(malt_x)
      sum_grad_sq=sum(grad_x^2)
      }
      v=rnorm(d)
      Delta=0
      for(j in 1:L){
        v=v-half*grad_x
        x=x+h*v
        grad_y=grad(x)
        Delta=Delta-half*sum(v*(grad_x+grad_y))
        v=v-half*grad_y
        v=eta*v+zeta*rnorm(d)
        grad_x=grad_y
      }
      if(i<=n/2){malt_x=x}else{
        Delta=small*(sum(grad_x^2)-sum_grad_sq)+U(x)-U_malt_x+Delta
        if(rexp(1)>Delta){malt_x=x}
        }
    }
    chain[1,]=malt_x
  }
  norm_draws=array(stats::rnorm(n*d*(L+1)),dim=c(n,d,L+1))
  expo_draws=stats::rexp(n)
  for(i in 1:(n-1)){
    x=malt_x
    grad_x=grad(malt_x)
    U_malt_x=U(malt_x)
    sum_grad_sq=sum(grad_x^2)
    v=norm_draws[i,,L+1]
    Delta=0
    for(j in 1:L){
      v=v-half*grad_x
      x=x+h*v
      grad_y=grad(x)
      Delta=Delta-half*sum(v*(grad_x+grad_y))
      v=v-half*grad_y
      v=eta*v+zeta*norm_draws[i,,j]
      grad_x=grad_y
    }
    Delta=small*(sum(grad_x^2)-sum_grad_sq)+U(x)-U_malt_x+Delta
    if(expo_draws[i]>Delta){malt_x=x}
    chain[i+1,]=malt_x
  }
  output=list(samples=chain, draw=chain[n,], acceptance=mean(chain[2:n,1]!=chain[1:(n-1),1]),g=g,h=h,L=L)
  class(output)=c("malt","mcmc")
  return(output)
  }
