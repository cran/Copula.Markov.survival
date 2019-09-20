#' Generate data from the Frank copula for serial dependence and the Frank copula for dependent censoring with the Weibull distribuions
#'
#' The data generation process is based on the Frank copula C_theta for serial dependence and the Frank copula tilde(C)_alpha for dependent censoring with the marginal disributions Weib(r, nu_1) and Weib(lambda, nu_2). Censoring percentage can be controled by constant c. This function is used when doing parametric bootstrap. The guide for using this function shall be explained by Huang (2019).
#'
#' @usage FrankFrank.Weibull.data(r, nu_1, theta, lambda, nu_2, alpha, N, c, l)
#'
#' @param r scale parameter for Weib(r, nu_1), r > 0
#' @param nu_1 shape parameter for Weib(r, nu_1), nu_1 > 0
#' @param theta copula parameter for C_theta, theta \eqn{\neq} 0
#' @param lambda scale parameter for Weib(lambda, nu_2), lambda > 0
#' @param nu_2 shape parameter for Weib(lambda, nu_2), nu_2 > 0
#' @param alpha copula parameter for tilde(C)_alpha, alpha \eqn{\neq} 0
#' @param N sample size
#' @param c parameter of Unif(0, c) for controlling censoring percentage
#' @param l length for data generation (default = 300)
#'
#' @return A list with the following elements:
#' \item{Subject}{a vector for numbers of subject}
#' \item{T_ij}{a vector for event times}
#' \item{delta_ij}{a vector for event indivator (=1 if recurrent; =0 if censoring)}
#' \item{T_i_star}{a vector for death times}
#' \item{delta_i_star}{a vector for death indivator (=1 if death; =0 if censoring)}
#'
#' @examples
#' Y = FrankFrank.Weibull.data(r = 1, nu_1 =0.5, theta = 2,
#'                             lambda = 0.45, nu_2 = 0.5, alpha = 2,
#'                             N = 100, c = 10, l = 300)
#'
#' @author Xinwei Huang
#'
#' @export

FrankFrank.Weibull.data = function(r, nu_1, theta, lambda, nu_2, alpha, N, c, l){

  U = matrix(NA, nrow = N, ncol = l)

  W = runif(n = N) #S_D(Di)
  A = runif(n = N) #random
  U[,1] = -1/alpha * log( ((exp(-alpha)-1)*A)/(exp(-alpha*W)-(exp(-alpha*W)-1)*A) + 1) #S_X(Xi1)
  for( j in c(2: l) ){
    A = runif(n = N) #random
    U[,j] = -1/theta * log( ((exp(-theta)-1)*A)/(exp(-theta*U[,j-1])-(exp(-theta*U[,j-1])-1)*A) + 1) #S_X(Xij)

  }

  Xij = ( -1/r * log(U) )^(1/nu_1) #Xij

  Ci = NA
  Di = NA
  T_i_star = NA
  delta_i_star = NA

  Ci = runif(N ,min = 0, max = c)
  Di = ( -1/lambda * log(W) )^(1/nu_2)
  T_i_star = pmin(Di,Ci)
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
