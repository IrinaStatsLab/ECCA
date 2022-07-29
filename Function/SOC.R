#----------------------------------
# Dongbang Yuan, Yunfeng Zhang
# Contact: yuandb09@gmail.com, yunfezhang@gmail.com
# Functions to update individual score matrices using splitting orthogonality constraints algorithm.
#------------------------------

#' Compute the updated Z1 and Z2 by SOC algorithm under orthogonal constraints when other parameters are fixed. 
#'
#' @param Z1 An n by r1 individual score matrix before update
#' @param Z2 An n by r2 individual score matrix before update
#' @param x1 First data matrix.
#' @param x2 Second data matrix.
#' @param family1 exponential family of x1.
#' @param family2 exponential family of x2.
#' @param U1 An n by r0 joint score matrix.
#' @param U2 An n by r0 joint score matrix.
#' @param V1 An p1 by r0 first joint loading matrix
#' @param V2 An p2 by r0 second joint loading matrix
#' @param A1 An p1 by r1 individual loading matrix
#' @param A2 An p2 by r2 individual loading matrix
#' @param mu1 The first mean vector
#' @param mu2 The second mean vector
#' @param max_iters Max iteration number for SOC main function.
#' @param newton_max_iters Max iteration number for damped Newton method.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param eps Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-4.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param alpha A scalar between 0 to 1. This parameter is for damped newton process and no need to change it.
#' @param gamma A positive number which is used in SOC update. Default is 1000. Enlarge this value if algorithm fails to convergence. If Algorithm is running too slow, pick a smaller number.
#' @param trials The parameter only used when the exponential family is binomial. It's the number of trials in binomial distribution. Default NULL.
#' 
#' @return A list containing the following output.
#' \itemize{
#'   \item Z1 - main output: updated first individual score matrix 
#'   \item Z2 - main output: updated second individual score matrix 
#'   \item B1 - B1 matrix in SOC algorithm 
#'   \item B2 - B2 matrix in SOC algorithm 
#'   \item P1 - P1 matrix in SOC algorithm 
#'   \item P2 - P2 matrix in SOC algorithm 
#'   \item obj - vector of negative log-likelihood recorded in the iteration steps
#'   \item h_vec - vector which monitors the constraints in the algorithm 
#'   \item criteria - a vector used to monitor convergence in SOC algorithm. The entry represents the maximum of dual and primal residuals at each iteration.
#' }
SOC <- function(Z1, Z2, x1, x2, family1, family2, U1, U2, V1, V2, A1, A2, mu1, mu2, max_iters = 1000, newton_max_iters = 100, 
                verbose = F, eps = 1e-4, learning_rate = 0.5, alpha = 0.25, gamma = 1000, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  r1 = ncol(A1)
  r2 = ncol(A2)
  n = nrow(x1)

  # record the fixed parts when computing log-likelihood. 
  offset1 = tcrossprod(rep(1, n), mu1) + tcrossprod(U1, V1)
  offset2 = tcrossprod(rep(1, n), mu2) + tcrossprod(U2, V2)
  
  # compute natural parameters
  theta1 = tcrossprod(Z1, A1) + offset1
  theta2 = tcrossprod(Z2, A2) + offset2
  
  # compute the projection matrix for U1 and U2 used in the while loop
  U = cbind(rep(1, n), U1, U2)
  projU = diag(n) - U %*% pinv(U)

  # negative log-likelihood
  # obj = - L1 - L2
  objective = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials)
  
  if (verbose){
    cat("The negative log-likelihood is ", objective, '\n')
  }
  
  # convergence criteria for SOC algorithm
  difference = 1
  # the vector including the convergence criteria at each step. 
  criteria = difference
  
  # The constraints. We only focus on constraints regarding to Z1 and Z2.
  h = max(Fnorm(crossprod(cbind(rep(1, n), U1, U2), cbind(Z1, Z2))), 
          Fnorm(crossprod(Z1, Z2)),
          Fnorm(crossprod(Z1, Z1) - diag(r1)),
          Fnorm(crossprod(Z2, Z2) - diag(r2)))
  h_vec = h
  # Initialize P and B for SOC algorithm 
  P1 = Z1
  P2 = Z2
  B1 = matrix(0, ncol = r1, nrow = n)
  B2 = matrix(0, ncol = r2, nrow = n)
  
  # Main while loop
  while (difference > eps | h > eps) {
    niter = niter + 1
    # Update Z1 and Z2
    if (family1 == "gaussian"){
      new_Z1 = Analyticalupdate(x1, mu1, U1, V1, B1, P1, A1, gamma)
   }
    else{
      new_Z1 = fitZ(Z1, offset1, A1, x1, family1, B1 - P1, gamma, 
                    max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    if (family2 == "gaussian"){
      new_Z2 = Analyticalupdate(x2, mu2, U2, V2, B2, P2, A2, gamma)
    }
    else{
      new_Z2 = fitZ(Z2, offset2, A2, x2, family2, B2 - P2, gamma, 
                    max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    # Update P1 & P2 by computing SVD
    svdP = svd(projU %*% (cbind(new_Z1, new_Z2) + cbind(B1, B2)))
    newP = tcrossprod(svdP$u, svdP$v)
    new_P1 = newP[, 1:r1, drop = F]
    new_P2 = newP[, (r1+1) : (r1+r2), drop = F]
    
    # Update B1 and B2
    new_B1 = B1 + new_Z1 - new_P1
    new_B2 = B2 + new_Z2 - new_P2
    
    # Compute updated likelihood
    theta1 = tcrossprod(new_Z1, A1) + offset1
    theta2 = tcrossprod(new_Z2, A2) + offset2
    new_obj = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials)
    objective = c(objective, new_obj)
    
    # Compute the constraints
    h = max(Fnorm(crossprod(cbind(rep(1, n), U1, U2), cbind(new_Z1, new_Z2))), 
            Fnorm(crossprod(new_Z1, new_Z2)),
            Fnorm(crossprod(new_Z1, new_Z1) - diag(r1)),
            Fnorm(crossprod(new_Z2, new_Z2) - diag(r2)))
    h_vec = c(h_vec,h)
    
    # Use ADMM criteria to monitor convergence. Calculate the maximum of dual and primal residual.
    dual_res = gamma*Fnorm(cbind(new_P1, new_P2) - cbind(P1, P2))
    primal_res = Fnorm(cbind(new_Z1, new_Z2) - cbind(new_P1, new_P2))
    difference = max(dual_res, primal_res)
    criteria = c(criteria, difference)
    
    if (verbose){
      cat("The objective function is ", new_obj, " The constraint is ", h, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    
    if (niter > max_iters){
      warning("SOC update of Z fails to converge.")
      break
    }
    # Prepare for the next iteration.
    Z1 = new_Z1
    Z2 = new_Z2
    P1 = new_P1
    P2 = new_P2
    B1 = new_B1
    B2 = new_B2
  }

  if(abs(h) > eps){
    warning("The constraint of SOC function is not satisfied.")
  }
  return(list(Z1 = Z1, Z2 = Z2, B1 = B1, B2 = B2, P1 = P1, P2 = P2, obj = objective, h_vec = h_vec, criteria = criteria))
}

#' Compute the updated Z1 and Z2 by SOC algorithm under orthogonal constraints when other parameters are fixed without main effect. 
#'
#' @param Z1 An n by r1 individual score matrix before update
#' @param Z2 An n by r2 individual score matrix before update
#' @param x1 First data matrix.
#' @param x2 Second data matrix.
#' @param family1 exponential family of x1.
#' @param family2 exponential family of x2.
#' @param U1 An n by r0 joint score matrix.
#' @param U2 An n by r0 joint score matrix.
#' @param V1 An p1 by r0 first joint loading matrix
#' @param V2 An p2 by r0 second joint loading matrix
#' @param A1 An p1 by r1 individual loading matrix
#' @param A2 An p2 by r2 individual loading matrix
#' @param max_iters Max iteration number for SOC main function.
#' @param newton_max_iters Max iteration number for damped Newton method.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param eps Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-4.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param alpha A scalar between 0 to 1. This parameter is for damped newton process and no need to change it.
#' @param gamma A positive number which is used in SOC update. Default is 1000. Enlarge this value if algorithm fails to convergence. If Algorithm is running too slow, pick a smaller number.
#' @param trials The parameter only used when the exponential family is binomial. It's the number of trials in binomial distribution. Default NULL.
#' 
#' @return A list containing the following output.
#' \itemize{
#'   \item Z1 - main output: updated first individual score matrix 
#'   \item Z2 - main output: updated second individual score matrix 
#'   \item B1 - B1 matrix in SOC algorithm 
#'   \item B2 - B2 matrix in SOC algorithm 
#'   \item P1 - P1 matrix in SOC algorithm 
#'   \item P2 - P2 matrix in SOC algorithm 
#'   \item obj - vector of negative log-likelihood recorded in the iteration steps
#'   \item h_vec - vector which monitors the constraints in the algorithm 
#'   \item criteria - a vector used to monitor convergence in SOC algorithm. The entry represents the maximum of dual and primal residuals at each iteration.
#' }
SOC_original <- function(Z1, Z2, x1, x2, family1, family2, U1, U2, V1, V2, A1, A2, max_iters = 1000, newton_max_iters = 100, 
                         verbose = F, eps = 1e-4, learning_rate = 0.5, alpha = 0.25, gamma = 1000, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  r1 = ncol(A1)
  r2 = ncol(A2)
  n = nrow(x1)
  
  # record the fixed parts when computing log-likelihood. 
  offset1 = tcrossprod(U1, V1)
  offset2 = tcrossprod(U2, V2)
  
  # compute natural parameters
  theta1 = tcrossprod(Z1, A1) + offset1
  theta2 = tcrossprod(Z2, A2) + offset2
  
  # compute the projection matrix for U1 and U2 used in the while loop
  U = cbind(U1, U2)
  projU = diag(n) - U %*% pinv(U)
  
  # negative log-likelihood
  objective = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials)
  
  if (verbose){
    cat("The negative log-likelihood is ", objective, '\n')
  }
  
  # The convergence criteria.
  difference = 1
  criteria = difference
  
  # The constraints. We only focus on constraints regarding to Z1 and Z2.
  h = max(Fnorm(crossprod(cbind(U1, U2), cbind(Z1, Z2))), 
          Fnorm(crossprod(Z1, Z2)),
          Fnorm(crossprod(Z1, Z1) - diag(r1)),
          Fnorm(crossprod(Z2, Z2) - diag(r2)))
  h_vec = h
  
  # augmented parameters for SOC algorithm 
  P1 = Z1
  P2 = Z2
  B1 = matrix(0, ncol = r1, nrow = n)
  B2 = matrix(0, ncol = r2, nrow = n)
  
  # Main while loop
  while (difference > eps | h > eps) {
    niter = niter + 1
    # Update Z1 & Z2
    if (family1 == "gaussian"){
      new_Z1 = Analyticalupdate(x1, mu = NULL, U1, V1, B1, P1, A1, gamma, main_effect = FALSE)
    }
    else{
      new_Z1 = fitZ(Z1, offset1, A1, x1, family1, B1 - P1, gamma, 
                    max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    if (family2 == "gaussian"){
      new_Z2 = Analyticalupdate(x2, mu = NULL, U2, V2, B2, P2, A2, gamma, main_effect = FALSE)
    }
    else{
      new_Z2 = fitZ(Z2, offset2, A2, x2, family2, B2 - P2, gamma, 
                    max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    # Update P1 & P2 by computing SVD
    svdP = svd(projU %*% (cbind(new_Z1, new_Z2) + cbind(B1, B2)))
    newP = tcrossprod(svdP$u, svdP$v)
    new_P1 = newP[, 1:r1, drop = F]
    new_P2 = newP[, (r1+1) : (r1+r2), drop = F]
    
    # Update B1 and B2
    new_B1 = B1 + new_Z1 - new_P1
    new_B2 = B2 + new_Z2 - new_P2
    
    # Compute updated likelihood
    theta1 = tcrossprod(new_Z1, A1) + offset1
    theta2 = tcrossprod(new_Z2, A2) + offset2
    new_obj = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials)
    objective = c(objective, new_obj)
    
    # Compute the constraints
    h = max(Fnorm(crossprod(cbind(U1, U2), cbind(new_Z1, new_Z2))), 
            Fnorm(crossprod(new_Z1, new_Z2)),
            Fnorm(crossprod(new_Z1, new_Z1) - diag(r1)),
            Fnorm(crossprod(new_Z2, new_Z2) - diag(r2)))
    h_vec = c(h_vec,h)
    
    # Use ADMM criteria to monitor convergence. Calculate the maximum of dual and primal residual.
    dual_res = gamma*Fnorm(cbind(new_P1, new_P2) - cbind(P1, P2))
    primal_res = Fnorm(cbind(new_Z1, new_Z2) - cbind(new_P1, new_P2))
    difference = max(dual_res, primal_res)
    criteria = c(criteria, difference)
    
    if (verbose){
      cat("The objective function is ", new_obj, " The constraint is ", h, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    if (niter > max_iters){
      warning("SOC update of Z fails to converge.")
      break
    }
    # Prepare for the next iteration.
    Z1 = new_Z1
    Z2 = new_Z2
    P1 = new_P1
    P2 = new_P2
    B1 = new_B1
    B2 = new_B2
  }
  if(abs(h) > eps){
    warning("The constraint of SOC function is not satisfied.")
  }
  return(list(Z1 = Z1, Z2 = Z2, B1 = B1, B2 = B2, P1 = P1, P2 = P2, obj = objective, h_vec = h_vec, criteria = criteria))
}


# Compute the updated Z1 or Z2 of the augmented Lagrangian form using damped Newton's method.
fitZ <- function(Z, offset, coefMat, x, family, BP, gamma, max_iters = 100, verbose = F, eps = 1e-6, learning_rate = 0.5, alpha = 0.25, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  p = ncol(x)
  n = nrow(x)
  r = ncol(coefMat)
  
  # compute natural parameters
  theta = tcrossprod(Z, coefMat) + offset
  
  # negative log-likelihood
  objective = -exp_fam_log_like(x, theta, family, trials) + gamma/2 * sum((Z + BP)^2)
  
  if (verbose){
    cat("The negative log-likelihood is ", objective, '\n')
  }
  
  # difference between estimators of each iteration
  difference = 1
  
  while (difference > eps) {
    niter = niter + 1
    # calculate coefficients
    theta = tcrossprod(Z, coefMat) + offset
    means = exp_fam_mean(theta, family, trials)
    weights = exp_fam_variance(theta, family, trials)
    eta = means - x
    new_Z = Z
    # Newton update for each row
    for(i in 1:n){
      gradient = crossprod(coefMat, eta[i, ]) + gamma * (Z[i, ] + BP[i, ])
      hessian  = crossprod(coefMat * weights[i, ], coefMat) + gamma * diag(r)
      # calculate the increment of Newton's method
      increment = c(solve(hessian, -gradient))
      # calculate damped newton parameter and its likelihood
      lambda = 1
      new_Z[i, ] = Z[i, ] + c(lambda * increment)
      theta = tcrossprod(new_Z, coefMat) + offset
      new_obj = -exp_fam_log_like(x, theta, family, trials) + gamma/2 * sum((new_Z + BP)^2)
      # Backtracking line search
      while (new_obj > objective[niter] + alpha * lambda * c(crossprod(gradient, increment))){
        lambda = lambda * learning_rate
        new_Z[i, ] = Z[i, ] + c(lambda * increment)
        theta = tcrossprod(new_Z, coefMat) + offset
        new_obj = -exp_fam_log_like(x, theta, family, trials) + gamma/2 * sum((new_Z + BP)^2)
      }
    }
    # compute updated likelihood
    difference = abs(new_obj - tail(objective, 1))
    objective = c(objective, new_obj)
    if (verbose){
      cat("The negative log-likelihood is ", new_obj, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    if (niter > max_iters){
      warning("Newton's method fails to converge.")
      break
    }
    Z = new_Z
  }
  return(Z)
}