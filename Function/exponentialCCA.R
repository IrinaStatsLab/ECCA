#----------------------------------
# Dongbang Yuan
# Contact: yuandb09@gmail.com
#------------------------------

require(pracma)
#' Apply ECCA algorithm to data from exponential family. Currently works for Gaussian and Binomial proportion data.
#' 
#' @param x1 An n by p1 input matrix.
#' @param x2 An n by p2 input matrix.
#' @param r0 joint rank. Column dimension for U1, U2.
#' @param r1 individual rank for x1. Column dimension for Z1. Caution: notation is different from the notation in the paper. Here r1 is the INDIVIDUAL rank.
#' @param r2 individual rank for x2. Column dimension for Z2. Caution: notation is different from the notation in the paper. Here r2 is the INDIVIDUAL rank.
#' @param family1,family2 Exponential family of x1 and x2. a choice of c("gaussian", "binomial", "poisson", "multinomial", "exponential")
#' @param main_effect Do you want to include the main effect \mu? Default is TRUE. 
#' @param trials The parameter only used in binomial case. It's the number of trials in binomial distribution. Default NULL.
#' @param rho The penalty parameter for rho/2 * ||U1 - U2||_F^2. Default 0. The choice of rho is under investigation.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param max_iters Max iteration number for the whole algorithm.
#' @param max_iters_newton Max iteration number for Newton's method.
#' @param max_iters_soc Max iteration number for SOC method.
#' @param mu1,mu2 The intercept vectors of length n. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param U1 An n by r0 joint score matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param U2 An n by r0 joint score matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param Z1 An n by r1 individual score matrix for x1. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param Z2 An n by r2 individual score matrix for x2. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param V1 An p1 by r0 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param V2 An p2 by r0 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param A1 An p1 by r1 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param A2 An p2 by r2 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param gamma A positive number which is used in SOC update. Default is 1000. Enlarge this value if algorithm fails to convergence. If Algorithm is running too slow, pick a smaller number.
#' @param alpha A scalar which is used in the backtracking line search of damped Newton's method. No need to change it. Default 0.25.
#' @param eps Threshold value used to determine the overall convergence of the optimization algorithm; the default value is 1e-04.
#' @param eps_newton Threshold value used to determine the convergence of the newton method; the default value is 1e-06.
#' @param eps_soc Threshold value used to determine the convergence of the SOC algorithm; the default value is 1e-04.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param M parameter for binomial distribution
#' @param testing Logical. If true, show the training process of inner functions.
#' @return A list containing the following output.
#' \itemize{
#'   \item mu1 - main output: first mean vector
#'   \item mu2 - main output: second mean vector 
#'   \item U1 - main output: first joint score matrix 
#'   \item V1 - main output: first joint loading matrix 
#'   \item U2 - main output: second joint score matrix 
#'   \item V2 - main output: second joint loading matrix
#'   \item Z1 - main output: first individual score matrix 
#'   \item A1 - main output: first individual loading matrix
#'   \item Z2 - main output: second individual score matrix 
#'   \item A2 - main output: second individual loading matrix
#'   \item iters - how many iterations are run
#'   \item loss_trace - vector of negative log-likelihood recorded in the iteration steps
#'   \item constraints - vector which monitors the constraints in the algorithm 
#'   \item gamma - the gamma used in SOC algorithm. 
#' }
#' @examples
#' @export
ecca_v2 <- function(x1, x2, r0, r1, r2, family1 = c("gaussian", "binomial", "poisson", "multinomial", "exponential"), 
                    family2 = c("gaussian", "binomial", "poisson", "multinomial", "exponential"), main_effect = TRUE, trials = NULL, rho = 0, 
                    learning_rate = 0.5, max_iters = 1000, max_iters_newton = 100, max_iters_soc = 1000, 
                    mu1, mu2, U1, U2, Z1, Z2, V1, V2, A1, A2, gamma = 1000, alpha = 0.25, 
                    eps = 1e-4, eps_newton = 1e-6, eps_soc = 1e-4, verbose = TRUE, M = 10, testing = F){
  if (!main_effect){
    result = ecca_original(x1 = x1, x2 = x2, r0 = r0, r1 = r1, r2 = r2, family1 = family1, family2 = family2, trials = trials, rho = rho, 
                           learning_rate = learning_rate, max_iters = max_iters, max_iters_newton = max_iters_newton, max_iters_soc = max_iters_soc, 
                           U1 = U1, U2 = U2, Z1 = Z1, Z2 = Z2, V1 = V1, V2 = V2, A1 = A1, A2 = A2, gamma = gamma, alpha = alpha,
                           eps = eps, eps_newton = eps_newton, eps_soc = eps_soc, verbose = verbose, M = M, testing = testing)
    return(result)
  }
  else{
    # check if input data agree with the given families
    family1 = match.arg(family1)
    check_family(x1, family1, trials)
    family2 = match.arg(family2)
    check_family(x2, family2, trials)  
    
    # Data pre-processing: record dimensions
    x1 = as.matrix(x1)
    x2 = as.matrix(x2)
    n = nrow(x1)
    if(n!= nrow(x2)){
      stop("Please check the dimension of your input data X.")
    }
    p1 = ncol(x1)
    p2 = ncol(x2)
    ones = rep(1, n)
    # Compute saturated natural parameter matrices
    saturated1 = saturated_para(x1, family1, M = M, trials = trials)
    centerdSat1 = scale(saturated1, scale = F)
    saturated2 = saturated_para(x2, family2, M = M, trials = trials)
    centerdSat2 = scale(saturated2, scale = F)
    # initialize mu if not given
    if (missing(mu1)) {
      mu1 = colMeans(saturated1)
    }
    
    if (missing(mu2)) {
      mu2 = colMeans(saturated2)
    }
    
    # initialize U, Z if not given
    if (missing(U1) || missing(U2) || missing(Z1) || missing(Z2)){
      theta1_est = as.matrix(svd(centerdSat1)$u[,1:(r0+r1),drop=F])
      theta2_est = as.matrix(svd(centerdSat2)$u[,1:(r0+r2),drop=F])
      result = angle_cal(theta1_est, theta2_est)
      U1 = result$principal_vector1[,1:r0,drop=F]
      U2 = result$principal_vector2[,1:r0,drop=F]
      U = cbind(U1,U2)
      # Initialize Z
      projUortho = diag(n) - U %*% pinv(U)
      Z1 = svd(projUortho %*% centerdSat1)$u[,1:r1,drop=F]
      projZ1ortho = diag(n) - tcrossprod(Z1)
      Z2 = svd(projZ1ortho %*% projUortho %*% centerdSat2)$u[,1:r2,drop=F]
    }
    
    # initialize loading matrices V, A if not given
    if (missing(V1) || missing(V2) || missing(A1) || missing(A2)){
      V1 = crossprod(centerdSat1,U1)
      V2 = crossprod(centerdSat2,U2)
      A1 = crossprod(centerdSat1,Z1)
      A2 = crossprod(centerdSat2,Z2)
    }
    
    # compute the objectives: likelihood + penalty
    theta1 = outer(ones, mu1) + tcrossprod(U1, V1) + tcrossprod(Z1, A1)
    theta2 = outer(ones, mu2) + tcrossprod(U2, V2) + tcrossprod(Z2, A2)
    loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
    
    # record the loss
    loss_trace = numeric(max_iters + 1)
    loss_trace[1] = loglike
    
    # Record the constraints
    constraint = cal_constraints(U1, U2, Z1, Z2)[[1]]
    constraints = numeric(max_iters + 1)
    constraints[1] = constraint
    
    # set counter and objective difference
    count = 0
    difference = 1
    
    if (verbose) {
      print(paste0("Number of run: ", count, ", negative log-likelihood: ", loss_trace[1], ". "))
    }
    
    while (abs(difference) > eps) {
      count = count + 1
      # record the fixed parts when compute log-likelihood. 
      loading1 = cbind(mu1,V1,A1)
      loading2 = cbind(mu2,V2,A2)
      
      # update mu1 & V1, A1
      if (verbose){
        print("Updating mu1, V1, A1")
      }
      
      # If the exponential family is Gaussian, we choose to use the analytical update for loading matrices.
      if (family1 == "gaussian"){
        S1 = cbind(ones, U1, Z1)
        temp1 = MainAnalyticalupdate(x1, S1, r0, r1)
        mu1 = temp1$mu
        V1 = temp1$V
        A1 = temp1$A
      }
      
      # If not Gaussian, use damped Newton's method
      else{
        temp1 = dampedNewton(x1, U1, Z1, loading1, family1, max_iters = max_iters_newton, verbose = testing, eps = eps_newton, learning_rate = learning_rate, alpha = alpha, trials = trials)
        mu1 = temp1[,1]
        V1 = temp1[,2:(1+r0), drop = F]
        A1 = temp1[,(2+r0):(1+r0+r1), drop = F]
      } 
      
      if (verbose){
        theta1 = outer(ones, mu1) + tcrossprod(U1, V1) + tcrossprod(Z1, A1)
        loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
        print(paste0("The negative log-likelihood after updating mu1, V1, A1 is ", loglike))
      }
      
      # update mu2 & V2, A2
      if (verbose){
        print("Updating mu2, V2, A2")
      }
      
      # If the exponential family is Gaussian, we choose to use the analytical update for loading matrices.
      if (family2 == "gaussian"){
        S2 = cbind(ones, U2, Z2)
        temp2 = MainAnalyticalupdate(x2, S2, r0, r2)
        mu2 = temp2$mu
        V2 = temp2$V
        A2 = temp2$A
      }
      
      # If not Gaussian, use damped Newton's method
      else{
        temp2 = dampedNewton(x2, U2, Z2, loading2, family2, max_iters = max_iters_newton, verbose = testing, eps = eps_newton, learning_rate = learning_rate, alpha = alpha, trials = trials)
        mu2 = temp2[,1]
        V2 = temp2[,2:(1+r0), drop = F]
        A2 = temp2[,(2+r0):(1+r0+r2), drop = F]
      }
      
      if (verbose){
        theta2 = outer(ones, mu2) + tcrossprod(U2, V2) + tcrossprod(Z2, A2)
        loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
        print(paste0("The negative log-likelihood after updating mu2, V2, A2 is ", loglike))
      }
      
      # update Z1 & Z2
      if (verbose){
        print("Updating Z1, Z2")
      }
      # update Z1 & Z2 by analytical update
      if (family1 == "gaussian" & family2 == "gaussian"){
        new_Z = AnalyticalSOC_Z(x1,x2,mu1,mu2,U1,U2,V1,V2,A1,A2)
        Z1 = new_Z$Z1
        Z2 = new_Z$Z2
      }
      # Usual SOC update
      else{    
        new_Z = SOC(Z1, Z2, x1, x2, family1, family2, U1, U2, V1, V2, A1, A2, mu1, mu2,
                    max_iters = max_iters_soc, newton_max_iters = max_iters_newton, verbose = testing, eps = eps_soc, 
                    learning_rate = learning_rate, alpha = alpha, gamma = gamma, trials = trials)
        Z1 = new_Z$Z1
        Z2 = new_Z$Z2
      }
      
      if (verbose){
        theta1 = outer(ones, mu1) + tcrossprod(U1, V1) + tcrossprod(Z1, A1)
        theta2 = outer(ones, mu2) + tcrossprod(U2, V2) + tcrossprod(Z2, A2)
        loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
        print(paste0("The negative log-likelihood after updating Z1, Z2 is ", loglike))
      }
      
      # update U1 & U2
      if (verbose){
        print("Updating U1, U2")
      }
      # Gaussian closed-form update.
      if (family1 == "gaussian" & family2 == "gaussian"){
        U1 = AnalyticalSOC_U(x1, mu1, Z1, Z2, V1, A1, first = TRUE)
        U2 = AnalyticalSOC_U(x2, mu2, Z1, Z2, V2, A2, first = FALSE)
      }
      # SOC update
      else{
        newU = optimizeU(U1, U2, x1, x2, family1, family2, Z1, Z2, V1, V2, A1, A2, mu1, mu2, rho = rho, 
                         max_iters = max_iters_soc, newton_max_iters = max_iters_newton, verbose = testing, eps = eps_soc,
                         learning_rate = learning_rate, alpha = alpha, gamma = gamma, trials = trials)
        U1 = newU$U1
        U2 = newU$U2
      }
      if (verbose){
        theta1 = outer(ones, mu1) + tcrossprod(U1, V1) + tcrossprod(Z1, A1)
        theta2 = outer(ones, mu2) + tcrossprod(U2, V2) + tcrossprod(Z2, A2)
        loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
        print(paste0("The negative log-likelihood after updating U1, U2 is ", loglike))
      }
      
      ## normalize U, V such that U satisfies the orthogonal constraint and the negative log-likelihood stays the same.
      if (verbose){
        print("Normalizing U1, U2")
      }
      svdUU  = svd(crossprod(U1,U2))
      U1 = U1 %*% svdUU$u
      U2 = U2 %*% svdUU$v
      V1 = V1 %*% svdUU$u
      V2 = V2 %*% svdUU$v
      
      # compute new log-likelihood
      theta1 = outer(ones, mu1) + tcrossprod(U1, V1) + tcrossprod(Z1, A1)
      theta2 = outer(ones, mu2) + tcrossprod(U2, V2) + tcrossprod(Z2, A2)
      loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
      
      if (verbose){
        print(paste0("The negative log-likelihood after normalizing U1, U2 is ", loglike))
        print(paste0("Number of run: ", count, ", negative log-likelihood: ", loglike, ". "))
      }
      
      # Record the objective values and the constraints. Prepare for the next iteration.
      loss_trace[count + 1] = loglike
      constraint = cal_constraints(U1, U2, Z1, Z2)[[1]]
      constraints[count + 1] = constraint
      
      difference = abs(loss_trace[count] - loss_trace[count + 1])
      
      if (count > max_iters){
        warning("Algorithm fails to converge.")
        break
      }
    }
    
    object <- list(mu1 = mu1,
                   mu2 = mu2,
                   U1 = U1,
                   V1 = V1,
                   U2 = U2,
                   V2 = V2,
                   Z1 = Z1,
                   A1 = A1,
                   Z2 = Z2,
                   A2 = A2,
                   iters = count,
                   loss_trace = loss_trace[1:(count + 1)],
                   constraints = constraints[1:(count + 1)],
                   gamma = gamma)
    return(object)
  }
}

#' Apply ECCA algorithm to data from exponential family without main effect. Currently works for Gaussian and Binomial proportion data.
#' 
#' @param x1 An n by p1 input matrix.
#' @param x2 An n by p2 input matrix.
#' @param r0 joint rank. Column dimension for U1, U2.
#' @param r1 individual rank for x1. Column dimension for Z1. Caution: notation is different from the notation in the paper. Here r1 is the INDIVIDUAL rank.
#' @param r2 individual rank for x2. Column dimension for Z2. Caution: notation is different from the notation in the paper. Here r2 is the INDIVIDUAL rank.
#' @param family1,family2 Exponential family of x1 and x2. a choice of c("gaussian", "binomial", "poisson", "multinomial", "exponential")
#' @param trials The parameter only used in binomial case. It's the number of trials in binomial distribution. Default NULL.
#' @param rho The penalty parameter for rho/2 * ||U1 - U2||_F^2. Default 0. The choice of rho is under investigation.
#' @param learning_rate A scalar between 0 to 1, which is used in the backtracking line search of damped Newton's method to shrink the step size. No need to change it. Default 0.5.
#' @param max_iters Max iteration number for the whole algorithm.
#' @param max_iters_newton Max iteration number for Newton's method.
#' @param max_iters_soc Max iteration number for SOC method.
#' @param U1 An n by r0 joint score matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param U2 An n by r0 joint score matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param Z1 An n by r1 individual score matrix for x1. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param Z2 An n by r2 individual score matrix for x2. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param V1 An p1 by r0 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param V2 An p2 by r0 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param A1 An p1 by r1 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param A2 An p2 by r2 loading matrix. If given, we use it to initialize ECCA. If not given, default initialization process is used.
#' @param gamma A positive number which is used in SOC update. Default is 1000. Enlarge this value if algorithm fails to convergence. If Algorithm is running too slow, pick a smaller number.
#' @param alpha A scalar which is used in the backtracking line search of damped Newton's method. No need to change it. Default 0.25.
#' @param eps Threshold value used to determine the overall convergence of the optimization algorithm; the default value is 1e-04.
#' @param eps_newton Threshold value used to determine the convergence of the newton method; the default value is 1e-06.
#' @param eps_soc Threshold value used to determine the convergence of the SOC algorithm; the default value is 1e-04.
#' @param verbose Logical. If False, the algorithm will stay silent. If True, it will print out fitting progress.
#' @param M parameter for binomial distribution
#' @param testing Logical. If true, show the training process of inner functions.
#' @return A list containing the following output.
#' \itemize{
#'   \item U1 - main output: first joint score matrix 
#'   \item V1 - main output: first joint loading matrix 
#'   \item U2 - main output: second joint score matrix 
#'   \item V2 - main output: second joint loading matrix
#'   \item Z1 - main output: first individual score matrix 
#'   \item A1 - main output: first individual loading matrix
#'   \item Z2 - main output: second individual score matrix 
#'   \item A2 - main output: second individual loading matrix
#'   \item iters - how many iterations are run
#'   \item loss_trace - vector of negative log-likelihood recorded in the iteration steps
#'   \item constraints - vector which monitors the constraints in the algorithm 
#'   \item gamma - the gamma used in SOC algorithm. 
#' }
#' @examples
ecca_original <- function(x1, x2, r0, r1, r2, family1 = c("gaussian", "binomial", "poisson", "multinomial", "exponential"), 
                          family2 = c("gaussian", "binomial", "poisson", "multinomial", "exponential"), trials = NULL, rho = 0, 
                          learning_rate = 0.5, max_iters = 1000, max_iters_newton = 100, max_iters_soc = 1000, 
                          U1, U2, Z1, Z2, V1, V2, A1, A2, gamma = 1000, alpha = 0.25, eps = 1e-4, 
                          eps_newton = 1e-6, eps_soc = 1e-4, verbose = TRUE, M = 10, testing = F){
  # check if input data agree with the given families
  family1 = match.arg(family1)
  check_family(x1, family1, trials)
  family2 = match.arg(family2)
  check_family(x2, family2, trials)  
  
  # Data pre-processing: record dimensions
  x1 = as.matrix(x1)
  x2 = as.matrix(x2)
  n = nrow(x1)
  if(n!= nrow(x2)){
    stop("Please check the dimension of your input data X.")
  }
  p1 = ncol(x1)
  p2 = ncol(x2)
  # Compute saturated natural parameter matrices
  saturated1 = saturated_para(x1, family1, M = M, trials = trials)
  saturated2 = saturated_para(x2, family2, M = M, trials = trials)
  # initialize U, Z if not given
  if (missing(U1) || missing(U2) || missing(Z1) || missing(Z2)){
    theta1_est = as.matrix(svd(saturated1)$u[,1:(r0+r1),drop=F])
    theta2_est = as.matrix(svd(saturated2)$u[,1:(r0+r2),drop=F])
    result = angle_cal(theta1_est, theta2_est)
    U1 = result$principal_vector1[,1:r0,drop=F]
    U2 = result$principal_vector2[,1:r0,drop=F]
    U = cbind(U1,U2)
    # Initialize Z
    projUortho = diag(n) - U %*% pinv(U)
    Z1 = svd(projUortho %*% saturated1)$u[,1:r1,drop=F]
    projZ1ortho = diag(n) - tcrossprod(Z1)
    Z2 = svd(projZ1ortho %*% projUortho %*% saturated2)$u[,1:r2,drop=F]
  }
  # initialize loading matrices V, A if not given
  if (missing(V1) || missing(V2) || missing(A1) || missing(A2)){
    V1 = crossprod(saturated1,U1)
    V2 = crossprod(saturated2,U2)
    A1 = crossprod(saturated1,Z1)
    A2 = crossprod(saturated2,Z2)
  }
  
  # compute the objectives: likelihood + penalty
  theta1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
  theta2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
  loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
  
  # record the loss
  loss_trace = numeric(max_iters + 1)
  loss_trace[1] = loglike
  
  # Record the constraints
  constraint = cal_constraints(U1, U2, Z1, Z2, main_effect = FALSE)[[1]]
  constraints = numeric(max_iters + 1)
  constraints[1] = constraint
  
  # set counter and objective difference
  count = 0
  difference = 1
  
  if (verbose) {
    print(paste0("Number of run: ", count, ", negative log-likelihood: ", loss_trace[1], ". "))
  }
  
  while (abs(difference) > eps) {
    count = count + 1
    # record the fixed parts when compute log-likelihood. 
    loading1 = cbind(V1,A1)
    loading2 = cbind(V2,A2)
    
    # update V1, A1
    if (verbose){
      print("Updating V1, A1")
    }
    
    # If the exponential family is Gaussian, we choose to use the analytical update for loading matrices.
    if (family1 == "gaussian"){
      S1 = cbind(U1, Z1)
      temp1 = MainAnalyticalupdate(x1, S1, r0, r1, main_effect = FALSE)
      V1 = temp1$V
      A1 = temp1$A
    }
    
    # If not Gaussian, use damped Newton's method
    else{
      temp1 = dampedNewton_original(x1, U1, Z1, loading1, family1, max_iters = max_iters_newton, verbose = testing, eps = eps_newton, learning_rate = learning_rate, alpha = alpha, trials = trials)
      V1 = temp1[,1:r0, drop = F]
      A1 = temp1[,(1+r0):(r0+r1), drop = F]
    } 
    
    if (verbose){
      theta1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
      loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
      print(paste0("The negative log-likelihood after updating V1, A1 is ", loglike))
    }
    
    # update V2, A2
    if (verbose){
      print("Updating V2, A2")
    }
    
    # If the exponential family is Gaussian, we choose to use the analytical update for loading matrices.
    if (family2 == "gaussian"){
      S2 = cbind(U2, Z2)
      temp2 = MainAnalyticalupdate(x2, S2, r0, r2, main_effect = FALSE)
      V2 = temp2$V
      A2 = temp2$A
    }
    
    # If not Gaussian, use damped Newton's method
    else{
      temp2 = dampedNewton_original(x2, U2, Z2, loading2, family2, max_iters = max_iters_newton, verbose = testing, eps = eps_newton, learning_rate = learning_rate, alpha = alpha, trials = trials)
      V2 = temp2[,1:r0, drop = F]
      A2 = temp2[,(1+r0):(r0+r2), drop = F]
    }
    
    if (verbose){
      theta2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
      loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
      print(paste0("The negative log-likelihood after updating V2, A2 is ", loglike))
    }
    
    # update Z1 & Z2
    if (verbose){
      print("Updating Z1, Z2")
    }
    # update Z1 & Z2 by analytical update
    if (family1 == "gaussian" & family2 == "gaussian"){
      new_Z = AnalyticalSOC_Z(x1,x2,mu1 = NULL,mu2 = NULL,U1,U2,V1,V2,A1,A2,main_effect=FALSE)
      Z1 = new_Z$Z1
      Z2 = new_Z$Z2
    }
    # Usual SOC update
    else{    
      new_Z = SOC_original(Z1, Z2, x1, x2, family1, family2, U1, U2, V1, V2, A1, A2,
                           max_iters = max_iters_soc, newton_max_iters = max_iters_newton, verbose = testing, eps = eps_soc,
                           learning_rate = learning_rate, alpha = alpha, gamma = gamma, trials = trials)
      Z1 = new_Z$Z1
      Z2 = new_Z$Z2
    }
    
    if (verbose){
      theta1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
      theta2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
      loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
      print(paste0("The negative log-likelihood after updating Z1, Z2 is ", loglike))
    }
    
    # update U1 & U2
    if (verbose){
      print("Updating U1, U2")
    }
    # If we prioritize the constraints, we will use the P's rather than U's as the new U. 
    # Otherwise, simply use the original SOC output
    if (family1 == "gaussian" & family2 == "gaussian"){
      U1 = AnalyticalSOC_U(x1, mu = NULL, Z1, Z2, V1, A1, first = TRUE, main_effect = FALSE)
      U2 = AnalyticalSOC_U(x2, mu = NULL, Z1, Z2, V2, A2, first = FALSE, main_effect = FALSE)
    }
    else{
      newU = optimizeU_original(U1, U2, x1, x2, family1, family2, Z1, Z2, V1, V2, A1, A2, rho = rho, 
                                max_iters = max_iters_soc, newton_max_iters = max_iters_newton, verbose = testing, eps = eps_soc, 
                                learning_rate = learning_rate, alpha = alpha, gamma = gamma, trials = trials)
      U1 = newU$U1
      U2 = newU$U2
    }
    if (verbose){
      theta1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
      theta2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
      loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
      print(paste0("The negative log-likelihood after updating U1, U2 is ", loglike))
    }
    
    ## normalize U, V such that U satisfies the orthogonal constraint and the negative log-likelihood stays the same.
    if (verbose){
      print("Normalizing U1, U2")
    }
    svdUU  = svd(crossprod(U1,U2))
    U1 = U1 %*% svdUU$u
    U2 = U2 %*% svdUU$v
    V1 = V1 %*% svdUU$u
    V2 = V2 %*% svdUU$v
    
    # compute new log-likelihood
    theta1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
    theta2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
    loglike <- -exp_fam_log_like(x1, theta1, family1, trials) - exp_fam_log_like(x2, theta2, family2, trials) + rho/2 * sum((U1 - U2)^2)
    
    if (verbose){
      print(paste0("The negative log-likelihood after normalizing U1, U2 is ", loglike))
      print(paste0("Number of run: ", count, ", negative log-likelihood: ", loglike, ". "))
    }
    
    # Record the objective values and the constraints. Prepare for the next iteration.
    loss_trace[count + 1] = loglike
    constraint = cal_constraints(U1, U2, Z1, Z2, main_effect = FALSE)[[1]]
    constraints[count + 1] = constraint
    
    difference = abs(loss_trace[count] - loss_trace[count + 1])
    
    if (count > max_iters){
      warning("Algorithm fails to converge.")
      break
    }
  }
  
  object <- list(U1 = U1,
                 V1 = V1,
                 U2 = U2,
                 V2 = V2,
                 Z1 = Z1,
                 A1 = A1,
                 Z2 = Z2,
                 A2 = A2,
                 iters = count,
                 loss_trace = loss_trace[1:(count + 1)],
                 constraints = constraints[1:(count + 1)],
                 gamma = gamma)
  return(object)
}