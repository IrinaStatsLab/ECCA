#----------------------------------
# Dongbang Yuan and Yunfeng Zhang
# Contact: yuandb09@gmail.com
# Functions to update joint score matrices using splitting orthogonality constraints algorithm.
#------------------------------

#' Compute the updated U1 and U2 by SOC algorithm, when other parameters are fixed.
#' 
#' @param U1 An n by r0 score matrix, before updates.
#' @param U2 An n by r0 score matrix, before updates.
#' @param x1 First data matrix.
#' @param x2 Second data matrix.
#' @param family1 exponential family of x1.
#' @param family2 exponential family of x2.
#' @param Z1 An n by r1 individual score matrix
#' @param Z2 An n by r2 individual score matrix
#' @param V1 An p1 by r0 joint loading matrix
#' @param V2 An p2 by r0 joint loading matrix
#' @param A1 An p1 by r1 individual loading matrix
#' @param A2 An p2 by r2 individual loading matrix
#' @param mu1 The first intercept vector
#' @param mu2 The second intercept vector
#' @param rho Penalty parameter for rho/2 * ||U1 - U2||_F^2, default 0. This input is still under investigation. Do not change it.
#' @param max_iters Max iteration number for SOC.
#' @param newton_max_iters Max iteration number for damped Newton's method.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param eps Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-4.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param alpha A scalar  between 0 to 1. This parameter is for damped newton process and no need to change it.
#' @param gamma A positive number which is used in SOC update. Default is 1000. Enlarge this value if algorithm fails to convergence. If Algorithm is running too slow, pick a smaller number.
#' @param trials The parameter only used in binomial case. It's the number of trials in binomial distribution. Default NULL.
#' 
#' @return A list containing the following output.
#' \itemize{
#'   \item U1 - main output: updated first joint score matrix 
#'   \item U2 - main output: updated second joint score matrix 
#'   \item B1 - B1 matrix in SOC algorithm 
#'   \item B2 - B2 matrix in SOC algorithm 
#'   \item P1 - P1 matrix in SOC algorithm 
#'   \item P2 - P2 matrix in SOC algorithm 
#'   \item obj - vector of negative log-likelihood recorded in the iteration steps
#'   \item h_vec - vector which monitors the constraints in the algorithm 
#'   \item criteria - a vector used to monitor convergence in SOC algorithm. The entry represents the maximum of dual and primal residuals at each iteration.
#' }
optimizeU <- function(U1, U2, x1, x2, family1, family2, Z1, Z2, V1, V2, A1, A2, mu1, mu2, rho = 0, max_iters = 1000, newton_max_iters = 100, 
                      verbose = F, eps = 1e-4, learning_rate = 0.5, alpha = 0.25, gamma = 1000, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  p1 = ncol(x1)
  p2 = ncol(x2)
  r0 = ncol(U1)
  n = nrow(x1)

  # record the fixed parts when computing log-likelihood. 
  offset1 = tcrossprod(rep(1, n), mu1) + tcrossprod(Z1, A1)
  offset2 = tcrossprod(rep(1, n), mu2) + tcrossprod(Z2, A2)
  
  # compute the projection matrix for Z1 and Z2
  Z = cbind(rep(1, n), Z1, Z2)
  projZ = diag(n) - Z %*% pinv(Z)

  # compute natural parameters
  theta1 = tcrossprod(U1, V1) + offset1
  theta2 = tcrossprod(U2, V2) + offset2
  
  # negative log-likelihood plus penalty
  # obj = - L1 - L2 + rho/2 * ||U1 - U2||_F^2
  objective = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
  criteria = 1
  
  if (verbose){
    cat("The negative log-likelihood is ", objective, '\n')
  }
  
  # difference between estimators of each iteration
  difference = 1
  
  # The constraints. We only focus on constraints regarding to U1 and U2.
  h = max(Fnorm(crossprod(cbind(rep(1, n),Z),cbind(U1,U2))),
          Fnorm(crossprod(U1, U1) - diag(r0)),
          Fnorm(crossprod(U2, U2) - diag(r0)))
  h_vec = h
  
  # augmented parameters for SOC algorithm 
  P1 = U1 
  P2 = U2
  B1 = matrix(0, ncol = ncol(U1), nrow = n)
  B2 = matrix(0, ncol = ncol(U2), nrow = n)
  while (difference > eps | h > eps){
    niter = niter + 1
    # Update U1 & U2
    if (family1 == "gaussian"){
      new_U1 = Analyticalupdate(x1, mu1, Z1, A1, B1, P1, V1, gamma)
    }
    else{
      new_U1 = convexFitU(U1, U2, offset1, V1, x1, family1, B1 - P1, gamma, rho = 0,
                          max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    if (family2 == "gaussian"){
      new_U2 = Analyticalupdate(x2, mu2, Z2, A2, B2, P2, V2, gamma)
    }
    else{
      new_U2 = convexFitU(U2, U1, offset2, V2, x2, family2, B2 - P2, gamma, rho = 0,
                          max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    ## Update P1 and P2. 
    svdP1 = svd(projZ %*% (new_U1 + B1))
    svdP2 = svd(projZ %*% (new_U2 + B2))
    new_P1 =  tcrossprod(svdP1$u, svdP1$v)
    new_P2 =  tcrossprod(svdP2$u, svdP2$v)
    
    # Update B1 and B2
    new_B1 = B1 + new_U1 - new_P1
    new_B2 = B2 + new_U2 - new_P2
    
    # compute updated likelihood and orthogonal constraint values
    theta1 = tcrossprod(new_U1, V1) + offset1
    theta2 = tcrossprod(new_U2, V2) + offset2
    new_obj = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((new_U1 - new_U2)^2)
    
    # The constraints. We only focus on constraints regarding to U1 and U2.
    h = max(Fnorm(crossprod(cbind(rep(1, n),Z),cbind(new_U1,new_U2))),
            Fnorm(crossprod(new_U1, new_U1) - diag(r0)),
            Fnorm(crossprod(new_U2, new_U2) - diag(r0)))
    
    h_vec = c(h_vec,h)
    
    # Use ADMM criteria to monitor convergence. Calculate the maximum of dual and primal residual.
    dual_res = gamma*Fnorm(cbind(new_P1, new_P2) - cbind(P1, P2))
    primal_res = Fnorm(cbind(new_U1, new_U2) - cbind(new_P1, new_P2))
    difference = max(dual_res, primal_res)

    # Record the objective value.
    objective = c(objective, new_obj)
    
    # Record the convergence criteria
    criteria = c(criteria, difference)
    
    if (verbose){
      cat("The objective function is ", new_obj, " The constraint is ", h, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    if (niter >= max_iters){
      warning("optU fails to converge.")
      break
    }
    # Prepare for the next iteration.
    U1 = new_U1
    U2 = new_U2
    P1 = new_P1
    P2 = new_P2
    B1 = new_B1
    B2 = new_B2
  }
  if(abs(h) > eps)
    warning("The constraint of optU function is not satisfied.")
  return(list(U1 = U1, U2 = U2, B1 = B1, B2 = B2, P1 = P1, P2 = P2, obj = objective, h_vec = h_vec, criteria = criteria))
}

#' Compute the updated U1 and U2 by SOC algorithm without main effect, when other parameters are fixed.
#' 
#' @param U1 An n by r0 score matrix, before updates.
#' @param U2 An n by r0 score matrix, before updates.
#' @param x1 First data matrix.
#' @param x2 Second data matrix.
#' @param family1 exponential family of x1.
#' @param family2 exponential family of x2.
#' @param Z1 An n by r1 individual score matrix
#' @param Z2 An n by r2 individual score matrix
#' @param V1 An p1 by r0 joint loading matrix
#' @param V2 An p2 by r0 joint loading matrix
#' @param A1 An p1 by r1 individual loading matrix
#' @param A2 An p2 by r2 individual loading matrix
#' @param rho Penalty parameter for rho/2 * ||U1 - U2||_F^2, default 0. This input is still under investigation. Do not change it.
#' @param max_iters Max iteration number for SOC.
#' @param newton_max_iters Max iteration number for damped Newton's method.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param eps Threshold value used to determine the convergence of the optimization algorithm; the default value is 1e-4.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param alpha A scalar  between 0 to 1. This parameter is for damped newton process and no need to change it.
#' @param gamma A positive number which is used in SOC update. Default is 1000. Enlarge this value if algorithm fails to convergence. If Algorithm is running too slow, pick a smaller number.
#' @param trials The parameter only used in binomial case. It's the number of trials in binomial distribution. Default NULL.
#' 
#' @return A list containing the following output.
#' \itemize{
#'   \item U1 - main output: updated first joint score matrix 
#'   \item U2 - main output: updated second joint score matrix 
#'   \item B1 - B1 matrix in SOC algorithm 
#'   \item B2 - B2 matrix in SOC algorithm 
#'   \item P1 - P1 matrix in SOC algorithm 
#'   \item P2 - P2 matrix in SOC algorithm 
#'   \item obj - vector of negative log-likelihood recorded in the iteration steps
#'   \item h_vec - vector which monitors the constraints in the algorithm 
#'   \item criteria - a vector used to monitor convergence in SOC algorithm. The entry represents the maximum of dual and primal residuals at each iteration.
#' }
optimizeU_original <- function(U1, U2, x1, x2, family1, family2, Z1, Z2, V1, V2, A1, A2, rho = 0, max_iters = 1000, newton_max_iters = 100, 
                               verbose = F, eps = 1e-4, learning_rate = 0.5, alpha = 0.25, gamma = 1000, trials = NULL)
{
  # record the dimensions and ranks
  niter = 0
  p1 = ncol(x1)
  p2 = ncol(x2)
  r0 = ncol(U1)
  n = nrow(x1)
  
  # record the fixed parts when computing log-likelihood. 
  offset1 = tcrossprod(Z1, A1)
  offset2 = tcrossprod(Z2, A2)
  
  # compute the projection matrix for Z1 and Z2
  Z = cbind(Z1, Z2)
  projZ = diag(n) - Z %*% pinv(Z)
  
  # compute natural parameters
  theta1 = tcrossprod(U1, V1) + offset1
  theta2 = tcrossprod(U2, V2) + offset2
  
  # negative log-likelihood plus penalty
  # obj = - L1 - L2 + rho/2 * ||U1 - U2||_F^2
  objective = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
  criteria = 1
  
  if (verbose){
    cat("The negative log-likelihood is ", objective, '\n')
  }
  
  # difference between estimators of each iteration
  difference = 1
  
  # The constraints. We only focus on constraints regarding to U1 and U2.
  h = max(Fnorm(crossprod(Z,cbind(U1,U2))), 
          Fnorm(crossprod(U1, U1) - diag(r0)),
          Fnorm(crossprod(U2, U2) - diag(r0)))
  h_vec = h
  
  # augmented parameters for SOC algorithm 
  P1 = U1 
  P2 = U2
  B1 = matrix(0, ncol = ncol(U1), nrow = n)
  B2 = matrix(0, ncol = ncol(U2), nrow = n)
  while (difference > eps | h > eps){
    niter = niter + 1
    # Update U1 & U2
    if (family1 == "gaussian"){
      new_U1 = Analyticalupdate(x1, mu1 = NULL, Z1, A1, B1, P1, V1, gamma, main_effect = FALSE)
    }
    else{
      new_U1 = convexFitU(U1, U2, offset1, V1, x1, family1, B1 - P1, gamma, rho = 0,
                          max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    if (family2 == "gaussian"){
      new_U2 = Analyticalupdate(x2, mu2 = NULL, Z2, A2, B2, P2, V2, gamma, main_effect = FALSE)
    }
    else{
      new_U2 = convexFitU(U2, U1, offset2, V2, x2, family2, B2 - P2, gamma, rho = 0,
                          max_iters = newton_max_iters, verbose = verbose, eps = eps, learning_rate = learning_rate, alpha = alpha, trials = trials)
    }
    ## Update P1 and P2. 
    svdP1 = svd(projZ %*% (new_U1 + B1))
    svdP2 = svd(projZ %*% (new_U2 + B2))
    new_P1 = tcrossprod(svdP1$u, svdP1$v)
    new_P2 = tcrossprod(svdP2$u, svdP2$v)
    
    # Update B1 and B2
    new_B1 = B1 + new_U1 - new_P1
    new_B2 = B2 + new_U2 - new_P2
    
    # compute updated likelihood and orthogonal constraint values
    theta1 = tcrossprod(new_U1, V1) + offset1
    theta2 = tcrossprod(new_U2, V2) + offset2
    new_obj = -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((new_U1 - new_U2)^2)
 
    # The constraints. We only focus on constraints regarding to U1 and U2.
    h = max(Fnorm(crossprod(Z,cbind(new_U1,new_U2))), 
            Fnorm(crossprod(new_U1, new_U1) - diag(r0)),
            Fnorm(crossprod(new_U2, new_U2) - diag(r0)))
    h_vec = c(h_vec,h)
    
    # Use ADMM criteria to monitor convergence. Calculate the maximum of dual and primal residual.
    dual_res = gamma*Fnorm(cbind(new_P1, new_P2) - cbind(P1, P2))
    primal_res = Fnorm(cbind(new_U1, new_U2) - cbind(new_P1, new_P2))
    difference = max(dual_res, primal_res)
    
    # Record the objective value.
    objective = c(objective, new_obj)
    
    # Record the convergence criteria
    criteria = c(criteria, difference)
    
    if (verbose){
      cat("The objective function is ", new_obj, " The constraint is ", h, '\n')
      print(paste0("Number of run: ", niter, ", difference: ", difference, ". "))
    }
    if (niter >= max_iters){
      warning("optU fails to converge.")
      break
    }
    # Prepare for the next iteration.
    U1 = new_U1
    U2 = new_U2
    P1 = new_P1
    P2 = new_P2
    B1 = new_B1
    B2 = new_B2
  }
  if(abs(h) > eps)
    warning("The constraint of optU function is not satisfied.")
  return(list(U1 = U1, U2 = U2, B1 = B1, B2 = B2, P1 = P1, P2 = P2, obj = objective, h_vec = h_vec, criteria = criteria))
}


# Compute the updated U1 OR U2 of the augmented Lagrangian form using damped Newton's method.
# The first argument is the part that will be updated and returned. 
# The second input U2 is only used for the 'unused' penalty term 'rho/2 * ||U1 - U2||_F^2'
# Arguments 'offset, coefMat, x, family, BP' are all corresponding to the first argument.
convexFitU <- function(U1, U2, offset, coefMat, x, family, BP, gamma, rho = 0, 
                       max_iters = 100, verbose = F, eps = 1e-6, learning_rate = 0.5, alpha = 0.25, trials = NULL)
{
  # record parameters
  niter = 0
  p = ncol(x)
  n = nrow(x)
  r = ncol(coefMat)
  
  # compute natural parameters
  theta = tcrossprod(U1, coefMat) + offset
  
  # objective values
  objective = -exp_fam_log_like(x, theta, family, trials) + gamma/2 * sum((U1 + BP)^2) + rho/2 * sum((U1 - U2)^2)
  
  if (verbose){
    cat("The negative log-likelihood of convexFit is ", objective, '\n')
  }
  
  # difference between objective values in each iteration
  difference = 1
  
  while (difference > eps) {
    niter = niter + 1
    # calculate coefficients
    theta = tcrossprod(U1, coefMat) + offset
    means = exp_fam_mean(theta, family, trials)
    weights = exp_fam_variance(theta, family, trials)
    eta = means - x
    new_U1 = U1
    # Newton update for each row
    for(i in 1:n){
      gradient = crossprod(coefMat, eta[i, ]) + gamma * (U1[i, ] + BP[i, ]) + rho * (U1[i, ] - U2[i, ])
      hessian  = crossprod(coefMat * weights[i, ], coefMat) + (gamma + rho) * diag(r) 
      # calculate the increment of Newton's method
      increment = c(solve(hessian, -gradient))
      # calculate damped newton parameter and corresponding likelihood
      lambda = 1
      new_U1[i, ] = U1[i, ] + c(lambda * increment)
      theta = tcrossprod(new_U1, coefMat) + offset
      new_obj = -exp_fam_log_like(x, theta, family, trials) + gamma/2 * sum((new_U1 + BP)^2) + rho/2 * sum((new_U1 - U2)^2)
      
      # choose step size
      while (new_obj > objective[niter] + alpha * lambda * c(crossprod(gradient, increment))){
        lambda = lambda * learning_rate
        new_U1[i, ] = U1[i, ] + c(lambda * increment)
        theta = tcrossprod(new_U1, coefMat) + offset
        new_obj = -exp_fam_log_like(x, theta, family, trials) + gamma/2 * sum((new_U1 + BP)^2) + rho/2 * sum((new_U1 - U2)^2)
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
    U1 = new_U1
  }
  return(U1)
}