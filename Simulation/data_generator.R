#----------------------------------
# Dongbang Yuan & Yunfeng Zhang
# Contact: yuandb09@gmail.com, yunfezhang@gmail.com
# Generate the data following ECCA model
#------------------------------

library(pracma)

# Generate random orthogonal matrix (n by p (n >= p))
GenOrthoMatrix <- function(n, p = n, tol = 1e-1){
  if (n < p){stop("n must be greater than p")}
  # initialize the rank to be 0.
  rank = 0
  while (rank < p){
    # Generate a random Gaussian matrix
    temp = matrix(rnorm(n^2), nrow = n, ncol = n)
    # SVD on the random matrix
    svd_temp = svd(temp)
    # Calculate the rank of the random matrix
    rank = sum(svd_temp$d > tol)
    # If the rank is equal to p, stops.
  }
  # The output is the full rank u matrix of the random Gaussian matrix
  result = svd_temp$u[,1:p]
  return(result)
}

# Center and Scale for each row in a matrix.
MatscaleRow <- function(X, center = TRUE, scale = TRUE){
  result = X
  if (center){
    Xmean = rowMeans(X)
    result = X - Xmean
  }
  n = dim(X)[1]
  p = dim(X)[2]
  sdvec = rep(0,n)
  if (scale){
    for (i in 1:n){
      tempstd = sd(X[i,])
      if (tempstd != 0 && !is.na(tempstd)){
        result[i,] = result[i,]/(tempstd * sqrt((p - 1) / p))
      }
      else{
        stop("The matrix either contains NA or there is no variation in some rows.")
        break
      }
    }
  }
  return(result)
}

# General function that does center and scale for a matrix, either row-wise or column-wise.
Matscale = function(X, center = TRUE, scale = TRUE, att = 'row'){
  if (att == 'row'){
    result = MatscaleRow(X, center = center, scale = scale)
  }
  if (att == 'col'){
    temp = MatscaleRow(t(X), center = center, scale = scale)
    result = t(temp)
  }
  return(result)
}

# Function to generate U.
Ugenerator <- function(n,r0,lambda, main_effect = TRUE){
  if (length(lambda)!=r0){
    stop("The length of lambda vector is not r0.")
  }
  covMat = rbind(cbind(diag(r0), diag(lambda, nrow = r0)), cbind(diag(lambda, nrow = r0), diag(r0)))
  
  # Generate joint scores
  G = matrix(rnorm(n^2), nrow = n, ncol = n) # Random matrix, n by n
  if (main_effect){
    G = Matscale(G, scale = FALSE, att = 'col') # center G, make sure our generated matrix is column centered.
  }
  U0 = svd(G)$u[,1:(2*r0), drop=F] # Get our random orthogonal U0, n by 2r0
  
  # Get the eigen-decomposition of covMat.
  eigen_result = eigen(covMat, symmetric = TRUE)
  
  # Slightly modify it to make sure our U is column centered, and crossprod(U) = covMat.
  U = U0 %*% diag(sqrt(eigen_result$values), nrow = 2*r0) %*% t(eigen_result$vectors)
  U1 = U[, 1:r0]
  U2 = U[, (r0+1):(2*r0)]
  return(list("U1" = U1, "U2" = U2))
}

# Function to generate natural parameters given U.
naturalParamGenerator = function(r0, r1, r2, n, U1, U2, p, family1, family2, mu_lb = 0.5, mu_ub = 1, sing_lb = 1, sing_ub = 1.2, scale_joint = 22, scale_ind1 = 15, scale_ind2 = 18, binomial_range = 3, trials = NULL, main_effect = TRUE){
  # Set singular values for joint and individual parts
  D0 = scale_joint*runif(r0,sing_lb,sing_ub)
  D1 = scale_ind1*runif(r1,sing_lb,sing_ub)
  D2 = scale_ind2*runif(r2,sing_lb,sing_ub)
  
  # Set dimensions
  p1 = p[1]
  p2 = p[2]
  
  # Generate intercept vectors
  if (main_effect){
    mu1 = runif(p1, mu_lb, mu_ub) * sample(c(1, -1), p1, T) # We put random +/- signal before the generated mean vectors.
    mu2 = runif(p2, mu_lb, mu_ub) * sample(c(1, -1), p2, T) 
    U = cbind(rep(1,n),U1, U2)
  }
  else{
    mu1 = NULL
    mu2 = NULL
    U = cbind(U1, U2)
  }
  # Generate individual scores
  Z = matrix(rnorm(n*(r1+r2),mean = 0, sd = 1), nrow = n, ncol = r1+r2)
  Z = svd((diag(n) - U %*% pinv(U)) %*% Z)$u
  Z1 = Z[, 1:r1,drop = F]
  Z2 = Z[, (r1 + 1):(r1 + r2), drop = F]
  # Generate loading matrices
  V1 = matrix(runif(p1 * r0, 1, 2) * sample(c(1, -1), p1 * r0, T), nrow = p1)
  V2 = matrix(runif(p2 * r0, 1, 2) * sample(c(1, -1), p2 * r0, T), nrow = p2)
  A1 = matrix(runif(p1 * r1, 1, 2) * sample(c(1, -1), p1 * r1, T), nrow = p1)
  A2 = matrix(runif(p2 * r2, 1, 2) * sample(c(1, -1), p2 * r2, T), nrow = p2)

  # Compute natural parameter matrices
  # If it's binomial case, we need to restrict our natural parameters (theta/n) to be roughly between (-3,3). theta/n = -3 means that p = 0.047; theta/n = 5 means that p = 0.99.
  if (family1 == 'binomial'){
    if (is.null(trials)){
      stop("You have specified the Binomial distribution. Please also tell us the number of trials.")
    }
    center1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
    # Control the scale of our centered natural parameter matrix to be in the range of (-3,3).
    scale_vec1 = rep(0,p1)
    for (i in 1:p1){
      # Get the half range and use it to scale the natural parameter.
      scale_vec1[i] = (max(center1[,i]) - min(center1[,i]))/2
      # Let the centered natural parameter to be within the range of (-3*trials, 3*trials)
      # The scaling won't change the joint or individual basis direction.
      center1[,i] = center1[,i]/scale_vec1[i] * binomial_range * trials
    }
    if (main_effect){
      theta1 = outer(rep(1, n), mu1) + center1
    }
    else{
      theta1 = center1
    }
    scaling_mat1 = diag(binomial_range * trials/scale_vec1, nrow = p1)
    V1 = scaling_mat1 %*% V1 # Re-scale V1 and A1 to make sure we have Theta1, Theta2 the model agrees.
    A1 = scaling_mat1 %*% A1
  }
  else{
    # Normalize the loading matrices
    V1 = svd(V1)$u
    A1 = svd(A1)$u
    V1 = V1 %*% diag(D0,nrow = r0, ncol = r0)
    A1 = A1 %*% diag(D1,nrow = r1, ncol = r1)
    if (main_effect){
      theta1 = outer(rep(1, n), mu1) + tcrossprod(U1, V1) + tcrossprod(Z1, A1)
    }
    else{
      theta1 = tcrossprod(U1, V1) + tcrossprod(Z1, A1)
    }
  }
  if (family2 == 'binomial'){
    if (is.null(trials)){
      stop("You have specified the Binomial distribution. Please also tell us the number of trials.")
    }
    center2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
    # Control the scale of our centered natural parameter matrix to be in the range of (-3,3).
    scale_vec2 = rep(0,p2)
    for (i in 1:p2){
      scale_vec2[i] = (max(center2[,i]) - min(center2[,i]))/2
      center2[,i] = center2[,i]/scale_vec2[i] * binomial_range * trials
    }
    if (main_effect){
      theta2 = outer(rep(1, n), mu2) + center2
    }
    else{
      theta2 = center2
    }
    scaling_mat2 = diag(binomial_range * trials/scale_vec2, nrow = p2)
    V2 = scaling_mat2 %*% V2
    A2 = scaling_mat2 %*% A2
  }
  else{
    V2 = svd(V2)$u
    A2 = svd(A2)$u
    V2 = V2 %*% diag(D0,nrow = r0, ncol = r0)
    A2 = A2 %*% diag(D2,nrow = r2, ncol = r2)
    if (main_effect){
      theta2 = outer(rep(1, n), mu2) + tcrossprod(U2, V2) + tcrossprod(Z2, A2)
    }
    else{
      theta2 = tcrossprod(U2, V2) + tcrossprod(Z2, A2)
    }
  }
  
  return (list(theta1 = theta1, 
               theta2 = theta2,
               mu1 = mu1,
               mu2 = mu2,
               V1 = V1,
               V2 = V2,
               Z1 = Z1,
               A1 = A1,
               Z2 = Z2,
               A2 = A2))
}

# This function generates datasets given natural parameters and exponential family.
# multinomial case has some issues.
# size is the number of trials in binomial cases.
# SNR is the ratio of Frobenious norm square between signals and noise.
dataSimulate = function(theta, family, size = 1e3, SNR = NULL, sd = NULL){
  n = nrow(theta)
  p = ncol(theta)
  if (is.null(SNR) & is.null(sd)){
    SNR = 1
  }
  # Calculate sd based on SNR, SNR priority is greater than sd.
  if (!is.null(SNR)){
    sd = Fnorm(theta)/sqrt((n*p*SNR))
  }
  if (family == 'gaussian'){
    # Normal case: natural parameters + Gaussian noise
    X = theta + matrix(rnorm(n * p, sd = sd), nrow = n)
  } 
  else if (family == "binomial") {
    # binomial case: natural parameters -> probability -> binomial variables
    probability = inv_logit_mat(theta/size)
    X = sapply(c(probability), function(x) rbinom(1, size, x)) /size
    X = matrix(X, nrow = n, ncol = p, byrow = FALSE)
  } 
  else if (family == "poisson") {
    # Poisson case
    lambdas = exp(theta)
    X = sapply(c(lambdas), function(x) rpois(1, lambda = x))
    X = matrix(X, nrow = n, ncol = p, byrow = FALSE)
  } 
  # else if (family == "multinomial") {
    # exp_theta = exp(theta)
    # mean_mat = sweep(exp_theta, 1, 1 + rowSums(exp_theta), "/")
  # } 
  else if (family == "exponential") {
    # Exponential case
    # Lambdas > 0 otherwise, rexp will give error.
    lambdas = -theta
    X = sapply(c(lambdas), function(x) rexp(1, rate = x))
    X = matrix(X, nrow = n, ncol = p, byrow = FALSE)
  }
  return(X)
}

# This function serves as a wrapper of the above functions. It generate datasets given ranks and exponential families.
dataGenerator = function(r0, r1, r2, n, p, lambda, family1, family2, size = 1e3, SNR = NULL, sd = NULL, main_effect = TRUE){
  # Generate natual parameters
  tempU = Ugenerator(n,r0,lambda, main_effect = main_effect)
  U1 = tempU$U1
  U2 = tempU$U2
  Params = naturalParamGenerator(r0, r1, r2, n, U1, U2, p, family1 = family1, family2 = family2, trials = size, main_effect = main_effect)
  # Generate datasets
  X1 = dataSimulate(Params[[1]], family1, size = size, SNR = SNR, sd = sd)
  X2 = dataSimulate(Params[[2]], family2, size = size, SNR = SNR, sd = sd)
  result = list(X1 = X1, X2 = X2, theta1 = Params$theta1, theta2 = Params$theta2, 
                U1 = U1, U2 = U2, Z1 = Params$Z1, Z2 = Params$Z2,
                mu1 = Params$mu1, mu2 = Params$mu2,
                V1 = Params$V1, V2 = Params$V2,
                A1 = Params$A1, A2 = Params$A2)
  return(result)
}
