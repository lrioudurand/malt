#' trajectory
#'
#' @description Draws a numerical Langevin trajectory with target density \deqn{\Pi(x)\propto e^{-U(x)}}{Pi(x)=exp(-U(x))/C} Given a potential function \eqn{U}{U} and its gradient evaluation, draws a trajectory and computes its numerical error. The trajectory drawn corresponds to the proposal in the sampling algorithm: Metropolis Adjusted Langevin Trajectories (Riou-Durand and Vogrinc 2022). Details available at: https://arxiv.org/abs/2202.13230.
#'
#' @param init Real vector. Initial values for the Langevin trajectory.
#' @param U A potential function to return the log-density of the distribution to be sampled from, up to an additive constant. It should input a real vector of the same length as init and output a scalar.
#' @param grad A function to return the gradient of the potential. It should input and output a real vector of the same length as init.
#' @param g Non-negative real number. The friction, a.k.a damping parameter. The choice g=0 boils down to Hamiltonian Monte Carlo.
#' @param h Positive real number. The time step.
#' @param L Positive integer. The number of steps per trajectory. The choice L=1 boils down to the Metropolis Adjusted Langevin Algorithm.
#'
#' @return Returns a list with the following objects:
#'   \item{path}{a matrix whose rows are the successive steps of the trajectory.}
#'   \item{draw}{a vector containing the output of the trajectory.}
#'   \item{num_error}{a vector containing the cumulative numerical errors along the path. The numerical error is measured by the energy difference incurred by the leapfrog integrator.}
#' @importFrom stats rnorm rexp
#' @export
#'
#' @examples
#' d=50
#' sigma=((d:1)/d)^(1/2)
#' init=rnorm(d)*sigma
#' U=function(x){sum(0.5*x^2/sigma^2)}
#' grad=function(x){x/sigma^2}
#' g=1.5
#' h=0.20
#' L=10
#' trajectory(init,U,grad,g,h,L)
trajectory=function(init,U,grad,g,h,L){
  d=length(init)
  if(d>=2){}else{stop("length(init) should be greater or equal than 2.")}
  if(g>=0){g=as.numeric(g)}else{warning("The friction/damping should be a non-negative real number, setting g=0 instead.");g=0}
  if(h>0){h=as.numeric(h)}else{stop("The time step should be a positive real number.")}
  if(L>=1){L=floor(L)}else{warning("The number of steps should be a positive integer, setting L=1 instead.");L=1}
  eta=exp(-g*h)
  zeta=sqrt(1-eta^2)
  path=matrix(NA,nrow=L+1,ncol=d)
  x=init
  path[1,]=x
  U_x=U(x)
  grad_x=grad(x)
  v=rnorm(d)
  Delta=rep(0,L+1)
  for(j in 1:L){
    kinetic=sum(v^2)
    v=v-(h/2)*grad_x
    x=x+h*v
    U_y=U(x)
    grad_y=grad(x)
    path[j+1,]=x
    v=v-(h/2)*grad_y
    kinetic_diff=(sum(v^2)-kinetic)/2
    Delta[j+1]=Delta[j]+U_y-U_x+kinetic_diff
    v=eta*v+zeta*stats::rnorm(d)
    grad_x=grad_y
    U_x=U_y
  }
  return(list(path=path,draw=path[L+1,],num_error=Delta))
}
