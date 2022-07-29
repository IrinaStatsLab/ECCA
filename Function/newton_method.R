#----------------------------------
# Dongbang Yuan
# Contact: yuandb09@gmail.com
# Modified from newton_method.R by Yunfeng Zhang.
#------------------------------

#' Damped Newton's method to update loading matrices. 
#' 
#' @param X Input matrix.
#' @param U The fixed joint score matrix U, of size n by r_0
#' @param Z The fixed individual score matrix Z, of size n by r_k (individual rank).
#' @param loading Initialized loading matrix cbind(mu,V,A).
#' @param family Exponential family of X.
#' @param max_iters Max iteration number. Default is 100.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param eps Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-06.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param alpha A scalar between 0 to 1. This parameter is for damped newton process and no need to change it.
#' @param trials The parameter is only used in binomial case. It's the number of trials in binomial distribution. Default NULL.
#' 
#' @return Updated column combined loading matrices cbind(mu, V, A).
dampedNewton <- function(X, U, Z, loading, family, max_iters = 100, verbose = F, eps = 1e-6, learning_rate = 0.5, alpha = 0.25, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  p = ncol(X)
  n = nrow(X)
  ones = matrix(rep(1,n), nrow = n)
  score = cbind(ones,U,Z)
  Theta = tcrossprod(score, loading)
  # negative log-likelihood
  llhood = -exp_fam_log_like(X, Theta, family, trials)
  
  if (verbose){
    cat("The negative log-likelihood is ", llhood, '\n')
  }
  
  # difference between estimators of each iteration
  difference = 1
  #Initialize new_loading matrix
  new_loading = loading
  while (difference > eps) {
    niter = niter + 1
    # calculate coefficients
    means = exp_fam_mean(Theta, family, trials)
    weights = exp_fam_variance(Theta, family, trials)
    eta = means - X
    # Update each row of loading matrix (mu,V,A)
    for(i in 1:p){
      gradient = crossprod(score, eta[, i])
      hessian  = crossprod(score * weights[, i], score)
      # calculate the increment of Newton's method
      increment = c(solve(hessian, -gradient))
      
      # calculate damped newton parameter and its likelihood
      lambda = 1
      new_loading[i, ] = loading[i, ] + c(lambda * increment)
      Theta = tcrossprod(score, new_loading)
      new_ll = -exp_fam_log_like(X, Theta, family, trials)
      # Backtracking line search to choose step size
      while (new_ll > llhood[niter] + alpha * lambda * c(crossprod(gradient, increment))){
        lambda = lambda * learning_rate
        new_loading[i, ] = loading[i, ] + c(lambda * increment)
        Theta = tcrossprod(score, new_loading)
        new_ll = -exp_fam_log_like(X, Theta, family, trials)
      }
    }
    # Check the difference of objective function.
    difference = abs(new_ll - tail(llhood, 1))
    llhood = c(llhood, new_ll)
    if (verbose){
      cat("The negative log-likelihood is ", new_ll, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    if (niter > max_iters){
      warning("Newton's method fails to converge.")
      break
    }
    loading = new_loading
  }
  if(llhood[1] - tail(llhood, 1) < -eps){
    warning("Newton's method likelihood increases.")}
  return(loading)
}

#' Damped Newton's method to update loading matrices without main effect. 
#' 
#' @param X Input matrix.
#' @param U The fixed joint score matrix U, of size n by r_0
#' @param Z The fixed individual score matrix Z, of size n by r_k (individual rank).
#' @param loading Initialized loading matrix cbind(V,A).
#' @param family Exponential family of X.
#' @param max_iters Max iteration number. Default is 100.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param eps Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-06.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param alpha A scalar between 0 to 1. This parameter is for damped newton process and no need to change it.
#' @param trials The parameter is only used in binomial case. It's the number of trials in binomial distribution. Default NULL.
#' 
#' @return Updated column combined loading matrices cbind(V, A)
dampedNewton_original <- function(X, U, Z, loading, family, max_iters = 100, verbose = F, eps = 1e-6, learning_rate = 0.5, alpha = 0.25, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  p = ncol(X)
  n = nrow(X)
  score = cbind(U,Z)
  Theta = tcrossprod(score, loading)
  # negative log-likelihood
  llhood = -exp_fam_log_like(X, Theta, family, trials)
  
  if (verbose){
    cat("The negative log-likelihood is ", llhood, '\n')
  }
  
  # difference between estimators of each iteration
  difference = 1
  #Initialize new_loading matrix
  new_loading = loading
  while (difference > eps) {
    niter = niter + 1
    # calculate coefficients
    means = exp_fam_mean(Theta, family, trials)
    weights = exp_fam_variance(Theta, family, trials)
    eta = means - X
    # Update each row of loading matrix (mu,V,A)
    for(i in 1:p){
      gradient = crossprod(score, eta[, i])
      hessian  = crossprod(score * weights[, i], score)
      # calculate the increment of Newton's method
      increment = c(solve(hessian, -gradient))
      
      # calculate damped newton parameter and its likelihood
      lambda = 1
      new_loading[i, ] = loading[i, ] + c(lambda * increment)
      Theta = tcrossprod(score, new_loading)
      new_ll = -exp_fam_log_like(X, Theta, family, trials)
      # Backtracking line search to choose step size
      while (new_ll > llhood[niter] + alpha * lambda * c(crossprod(gradient, increment))){
        lambda = lambda * learning_rate
        new_loading[i, ] = loading[i, ] + c(lambda * increment)
        Theta = tcrossprod(score, new_loading)
        new_ll = -exp_fam_log_like(X, Theta, family, trials)
      }
    }
    # Check the difference of objective function.
    difference = abs(new_ll - tail(llhood, 1))
    llhood = c(llhood, new_ll)
    if (verbose){
      cat("The negative log-likelihood is ", new_ll, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    if (niter > max_iters){
      warning("Newton's method fails to converge.")
      break
    }
    loading = new_loading
  }
  
  if(llhood[1] - tail(llhood, 1) < -eps){
    warning("Newton's method likelihood increases.")}

  return(loading)
}