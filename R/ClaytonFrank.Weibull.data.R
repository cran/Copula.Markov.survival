#' Generate data from the Clayton copula for serial dependence and the Frank copula for dependent censoring with the Weibull distributions
#'
#' The data generation process is based on the Clayton copula C_theta for serial dependence and the Frank copula tilde(C)_alpha for dependent censoring with the marginal distributions Weib(scale1, shape1) and Weib(scale2, shape2). Censoring percentage can be controlled by constant b. This function is used when doing parametric bootstrap. The guide for using this function shall be explained by Huang (2019), and Huang, Wang and Emura (2020).
#'
#' @usage ClaytonFrank.Weibull.data(N, scale1, shape1, theta, scale2, shape2, alpha, b, l)
#'
#' @param N sample size
#' @param scale1 scale parameter for Weib(scale1, shape1), scale1 > 0
#' @param shape1 shape parameter for Weib(scale1, shape1), shape1 > 0
#' @param theta copula parameter for C_theta, theta > 0
#' @param scale2 scale parameter for Weib(scale2, shape2), scale2 > 0
#' @param shape2 shape parameter for Weib(scale2, shape2), shape2 > 0
#' @param alpha copula parameter for tilde(C)_alpha, alpha \eqn{\neq} 0
#' @param b parameter of Unif(0, b) for controlling censoring percentage
#' @param l length for data generation (default = 300)
#'
#' @return A list with the following elements:
#' \item{Subject}{a vector for numbers of subject}
#' \item{T_ij}{a vector for event times}
#' \item{delta_ij}{a vector for event indicator (=1 if recurrent; =0 if censoring)}
#' \item{T_i_star}{a vector for death times}
#' \item{delta_i_star}{a vector for death indicator (=1 if death; =0 if censoring)}
#'
#' @examples
#' Y = ClaytonFrank.Weibull.data(N = 100, scale1 = 1, shape1 =0.5, theta = 2,
#'                               scale2 = 0.45, shape2 = 0.5, alpha = 2, b = 10, l = 300)
#'
#' @author Xinwei Huang
#' @references  Huang XW, Wang W, Emura T (2020) A copula-based Markov chain model for serially dependent event times with a dependent terminal event. Japanese Journal of Statistics & Data Science. Accepted.
#'
#' @export


ClaytonFrank.Weibull.data = function(N, scale1, shape1, theta, scale2, shape2, alpha, b, l){

  U = matrix(NA, nrow = N, ncol = l)

  W = runif(n = N) #S_D(Di)
  A = runif(n = N) #random
  U[,1] = -1/alpha * log( ((exp(-alpha)-1)*A)/(exp(-alpha*W)-(exp(-alpha*W)-1)*A) + 1) #S_X(Xi1)
  for( j in c(2: l) ){

    U[,j] = ((runif(n = N)^(-theta/(theta+1)) - 1) * U[,j-1]^(-theta) + 1)^(-1/theta) #S_X(Xij)

  }

  Xij = ( -1/scale1 * log(U) )^(1/shape1) #Xij

  Bi = NA
  Di = NA
  T_i_star = NA
  delta_i_star = NA

  Bi = runif(N ,min = 0, max = b)
  Di = ( -1/scale2 * log(W) )^(1/shape2)
  T_i_star = pmin(Di,Bi)
  delta_i_star = as.numeric(T_i_star == Di)

  Subject = T_ij = delta_ij = NULL

  for( i in c(1:N) ){

    cum.Xij = cumsum(Xij[i,])
    ni = min(which(cum.Xij > T_i_star[i]))

    if(ni == 1){
      Subject = c(Subject, i)
      T_ij = c(T_ij, T_i_star[i])
      delta_ij = c(delta_ij, 0)
    }else{
      Subject = c( Subject, rep(i,ni) )
      T_ij = c( T_ij, c(Xij[i,1:(ni-1)], T_i_star[i] - cum.Xij[ni-1]) )
      delta_ij = c( delta_ij, c(rep(1,ni-1),0) )
    }

  }

  return(list(Subject = Subject, T_ij = T_ij, delta_ij = delta_ij, T_i_star = T_i_star, delta_i_star = delta_i_star))

}
