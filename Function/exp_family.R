#----------------------------------
# Dongbang Yuan
# Contact: yuandb09@gmail.com
# Part of this file is a modification of this page https://github.com/andland/generalizedPCA/blob/master/R/exponential_family_functions.R
# This file includes the preliminary functions regarding to processing exponential family.
# This file also includes several functions for convenience.
#------------------------------

require(pracma)
# compute inverse logit function
inv_logit_mat <- function(x, min = 0, max = 1) {
  p <- exp(x)/(1 + exp(x))
  p <- ifelse(is.na(p) & !is.na(x), 1, p)
  p * (max - min) + min
}

# Check if the data fits the provided family
# family could be c('multinomial', 'binomial', 'poisson', 'exponential', 'gaussian') (no need to check gaussian).
# If violated, the reason will be shown.
# multinomial, poisson, exponential distributions are still under study. (As of 06/24/2022)
check_family <- function(x, family, trials = NULL) {
  distinct_vals = unique(c(x[!is.na(x)]))
  # Given family, check values (check binomial and multinomial)
  if ((family == "multinomial") & any(distinct_vals < 0 | distinct_vals > 1)) {
    stop(paste(family, "family with data outside [0, 1]"))
  }
  if (family == "binomial") {
    if (is.null(trials)){
      stop("You have specified the Binomial distribution. Please also tell us the number of trials.")
    }
    if (any(distinct_vals < 0 | distinct_vals/trials > 1)){
      stop(paste(family, "family with data outside [0, 1]"))
    }
  }
  # check poisson
  if (family == "poisson") {
    if (any(range(distinct_vals) < 0)) 
      stop("Negative data found for poisson family")
    if (any(distinct_vals%%1 != 0)) 
      warning("Non-integer data found for poisson family")
  }
  
  # Check exponential distribution
  if (family == "exponential") {
    if (any(range(distinct_vals) < 0)) 
      stop("Negative data found for exponential family")
  }
  
  if (family == "multinomial") {
    if (dim(x)[2] == 1){
      stop("The number of columns is 1. We cannot use the multinomial model.")
    }
    if (any(rowSums(x) != 1)) 
      stop("There is a row that do not sum up to 1. Please pre-process your data.")
    if (any(x[,ncol(x)] == 0)){
      stop("The last column contains 0 element. Please rearrange your data.")
    }
  }
  
  # Given data, check exponential family
  if (all(distinct_vals %in% c(0, 1))) {
    if (!(family %in% c("binomial", "multinomial"))) {
      message("All entries are binary. Do you want to use binomial?")
    }
  } else if (all(distinct_vals >= 0 & distinct_vals%%1 == 0)) {
    if (family != "poisson") {
      message("All entries are counts. Do you want to use poisson?")
    }
  }
}

# Given data and distribution, compute saturated natural parameters
# M is a user-given parameter that controls the scale of the natural parameter of binomial cases.
# If we observe x = 0, we get p = 0, which means eta = +Inf, which is tough to deal with when computing.
# E.g. if M = 10, then eta = -10m when we observe x = 0 with m trials.
# For Poisson distribution, if we observe x = 0, we set eta = -|M|.
# multinomial, poisson, exponential distributions are still under study. (As of 06/24/2022)
saturated_para <- function(x, family, M = 10, trials = NULL) {
  # Gaussian distribution
  if (family == "gaussian") {
    eta = x
  } 
  # Binomial distribution
  if (family == "binomial") {
    if (is.null(trials)){
      warning("You have specified the Binomial distribution. The result is based on only one trial. If you know the number of trials, please specify.")
    }
    # Initialize eta with |M| * trials or -|M| * trials.
    # Each entry in the data matrix comes from a Binomial distribution with entry of x/trials in (0,1) . 
    eta = abs(M) * (2 * x - 1) * trials
    non_binary = (x != 0 & x != 1 & !is.na(x))
    if (sum(non_binary) > 0) {
      # Calculate natural parameters for non-binary cases
      logitvals = log(x) - log(1 - x)
      eta[non_binary] = trials * logitvals[non_binary]
    }
  }
  # Multinomial distribution
  # Caution! If the input x is of size n by p
  # The output natural parameter matrix will have size n by p-1.
  # Because we treat each row as an outcome from a Multinomial distribution.
  # Refer to section 3.3 in the draft.
  if (family == "multinomial") {
    eta = abs(M) * (2 * x - 1)
    non_binary = (x != 0 & x != 1 & !is.na(x))
    if (sum(non_binary) > 0) {
      last_col = x[,ncol(x)]
      logitvals = log(x) - log(last_col)
      eta[non_binary] = logitvals[non_binary]
      eta = eta[,1:ncol(eta)-1]
    }
  }
  # Poisson distribution
  if (family == "poisson") {
    eta = log(x)
    eta[x==0] = -abs(M)
  }
  # Exponential distribution
  if (family == "exponential") {
    eta = -1/x
  }
  eta[is.na(x)] <- 0
  return(eta)
}

# compute the expectation of the random variable, which is mu=b'(theta)
# If the family is binomial, we also need the number of trials as an input. 
# multinomial, poisson, exponential distributions are still under study. (As of 06/24/2022)
exp_fam_mean <- function(theta, family, trials = NULL) {
  if (family == "gaussian") {
    mean_mat = theta
  } 
  else if (family == "binomial") {
    if (is.null(trials)){
      stop("You have specified the Binomial distribution. Please also tell us the number of trials.")
    }
    mean_mat = inv_logit_mat(theta/trials)
  } 
  else if (family == "poisson") {
    mean_mat = exp(theta)
  } 
  else if (family == "multinomial") {
    # If the family is multinomial, we will have the size of theta to be n by (p-1)
    exp_theta = exp(theta)
    mean_mat = exp_theta/(1 + rowSums(exp_theta))
  }
  else if (family == "exponential") {
    mean_mat = -1/theta
  }
  return(mean_mat)
}

# compute the variance, which is b''(theta)
# If the family is binomial, we also need the number of trials as an input.
# multinomial, poisson, exponential distributions are still under study. (As of 06/24/2022)
exp_fam_variance <- function(theta, family, trials = NULL) {
  if (family == "gaussian") {
    # The variance is 1 everywhere based on our assumptions.
    var_mat = matrix(1, nrow(theta), ncol(theta))
  }
  else if (family == "binomial") {
    if (is.null(trials)){
      stop("You have specified the Binomial distribution. Please also tell us the number of trials.")
    }
    mean_mat = exp_fam_mean(theta, family, trials)
    var_mat = 1/trials * mean_mat * (1 - mean_mat)
  }
  else if (family == "multinomial") {
    # If the family is multinomial, we will have the size of theta to be n by (p-1)
    mean_mat = exp_fam_mean(theta, family)
    var_mat = mean_mat * (1 - mean_mat)
  }
  else if (family == "poisson") {
    # Same as mean.
    var_mat = exp(theta)
  }
  else if (family == "exponential") {
    var_mat = 1/(theta^2)
  }
  return(var_mat)
}

#' Compute the log-likelihoods of given datasets.
#'
#' @param x the input dataset
#' @param theta natural parameters
#' @param family the exponential family of the dataset
#' 
#' @return the sum of the log-likelihoods of each observation.
#' @export
# multinomial, poisson, exponential distributions are still under study. (As of 06/24/2022)
exp_fam_log_like <- function(x, theta, family, trials = NULL) {
  if (family == "gaussian") {
    return(-0.5 * sum((x - theta)^2 , na.rm = TRUE))
  } 
  else if (family == "binomial") {
    if (is.null(trials)){
      stop("You have specified the Binomial distribution. Please also tell us the number of trials.")
    }
    # If there are large values in theta/trials, exp(theta/trials) will exceed the computation limit. 
    # In this case, log(1 + exp(theta/trials)) roughly equals to theta/trials 
    return(sum(x * theta - trials * log(1 + (theta/trials < 20) * exp(theta/trials)) - theta * (theta/trials >= 20), na.rm = TRUE))
  }
  else if (family == "poisson") {
    return(sum(x * theta - exp(theta) - lfactorial(x), na.rm = TRUE))
  } 
  else if (family == "multinomial") {
    if (any(rowSums(x) != 1)) 
      stop("There is a row that do not sum up to 1. Please check your data!")
    n = dim(x)[1]
    # Add back the last natural parameter column.
    theta = cbind(theta,rep(0,n))
    # Calculate the removed p_k. p_k is the probability when choosing the last category, for each row.
    last_p_vec = 1/rowSums(exp(theta))
    return(sum(x * theta - lfactorial(x), na.rm = TRUE) + sum(log(last_p_vec), na.rm = TRUE))
  }
  else if (family == "exponential") {
    return(sum(x * theta + log(-theta), na.rm = TRUE))
  }
}

# This function uses the adjustments mentioned in Dr. Longnecker's book to modify 0/1s in the binomial cases.
# If there are 0/1 in the binomial dataset, please run this function first to prevent malfunctions of ECCA.
# To be specific, 0 will be adjusted to 3/8/(n + 3/4)
# 1 will be adjusted to (n + 3/8)/(n + 3/4)
# n is the number of trials.
modifyBinomial <- function(X,trials){
  for (i in 1:dim(X)[1]){
    for (j in 1:dim(X)[2]){
      if (X[i,j] == 1){
        X[i,j] = (trials + 3/8)/(trials + 3/4)
      }
      if (X[i,j] == 0){
        X[i,j] = 3/8/(trials + 3/4)
      }
    }
  }
  return(X)
}

# Check if given joint and individual scoring matrices satisfy the constraints.
# Return Boolean variable
sanitycheck <- function(U1, U2, Z1, Z2, tol = 1e-6){
  n = dim(U1)[1]
  # Check whether U and Z are orthogonal
  if (sum(crossprod(U1,Z1)^2) > tol){
    warning("U1,Z1 failed.")
    return(0)
  }
  if (sum(crossprod(U1,Z2)^2) > tol){
    warning("U1,Z2 failed.")
    return(0)
  }
  if (sum(crossprod(U2,Z1)^2) > tol){
    warning("U2,Z1 failed.")
    return(0)
  }
  if (sum(crossprod(U2,Z2)^2) > tol){
    warning("U2,Z2 failed.")
    return(0)
  }
  # Check whether Z1 and Z2 are orthogonal
  if (sum(crossprod(Z1,Z2)^2) > tol){
    warning("Z1,Z2 failed.")
    return(0)
  }
  # Check whether U1, U2, Z1 and Z2 are column centered
  if (sum(colMeans(matrix(U1))^2) > tol){
    warning("Mean of U1 failed.")
    return(0)
  }
  if (sum(colMeans(matrix(U2))^2) > tol){
    warning("Mean of U2 failed.")
    return(0)
  }
  if (sum(colMeans(matrix(Z1))^2) > tol){
    warning("Mean of Z1 failed.")
    return(0)
  }
  if (sum(colMeans(matrix(Z2))^2) > tol){
    warning("Mean of Z2 failed.")
    return(0)
  }
  # Check whether the product of U1 and U2 is diagonal
  temp = crossprod(U1,U2)
  for (i in 1:dim(temp)[1]){
    for (j in 1:dim(temp)[2]){
      if (i!=j & abs(temp[i,j])>tol){
        warning("Cross product of U1 and U2 is not diagnonal.")
        return(0)
      }
    }
  }
  print("Sanity check passed!")
  return(1)
}

# Calculate Frobenius norm of a matrix
Fnorm <- function(X){
  X = as.matrix(X)
  result = sqrt(sum(X^2))
  return(result)
}

# A function to quantify the extent of the constraints satisfied by U1, U2, Z1 and Z2.
# Also return a message determining which constraint is controlling the extent
cal_constraints <- function(U1, U2, Z1, Z2, main_effect = TRUE){
  res = -Inf
  n = dim(U1)[1]
  # Check whether U and Z are orthogonal
  res = max(res, max(abs(crossprod(U1,Z1))))
  message = "The quantified constraints are provided by U1 and Z1."
  if (max(abs(crossprod(U1,Z2))) > res){
    message = "The quantified constraints are provided by U1 and Z2."
    res = max(abs(crossprod(U1,Z2)))
  }
  if (max(abs(crossprod(U2,Z1))) > res){
    message = "The quantified constraints are provided by U2 and Z1."
    res = max(abs(crossprod(U2,Z1)))
  }
  if (max(abs(crossprod(U2,Z2))) > res){
    message = "The quantified constraints are provided by U2 and Z2."
    res = max(abs(crossprod(U2,Z2)))
  }
  # Check whether Z1 and Z2 are orthogonal
  if (max(abs(crossprod(Z1,Z2))) > res){
    message = "The quantified constraints are provided by Z1 and Z2."
    res = max(abs(crossprod(Z1,Z2)))
  }
 
  # Check whether U1, U2, Z1 and Z2 are column centered
  if (main_effect){
    if (max(abs(colSums(U1))) > res){
      message = "The quantified constraints are provided by the orthognality between 1 and U1."
      res = max(abs(colSums(U1)))
    }
    if (max(abs(colSums(U2))) > res){
      message = "The quantified constraints are provided by the orthognality between 1 and U2."
      res = max(abs(colSums(U2)))
    }
    if (max(abs(colSums(Z1))) > res){
      message = "The quantified constraints are provided by the orthognality between 1 and Z1."
      res = max(abs(colSums(Z1)))
    }
    if (max(abs(colSums(Z2))) > res){
      message = "The quantified constraints are provided by the orthognality between 1 and Z2."
      res = max(abs(colSums(Z2)))
    }
  }
  r0 = dim(U1)[2]
  r1 = dim(Z1)[2]
  r2 = dim(Z2)[2]
  # Check whether U1, U2, Z1 and Z2 are orthogonal
  if (max(abs(crossprod(U1) - diag(r0))) > res){
    message = "The quantified constraints are provided by the orthogonality of U1."
    res = max(abs(crossprod(U1) - diag(r0)))
  }
  if (max(abs(crossprod(U2) - diag(r0))) > res){
    message = "The quantified constraints are provided by the orthogonality of U2."
    res = max(abs(crossprod(U2) - diag(r0)))
  }
  if (max(abs(crossprod(Z1) - diag(r1))) > res){
    message = "The quantified constraints are provided by the orthogonality of Z1."
    res = max(abs(crossprod(Z1) - diag(r1)))
  }
  
  if (max(abs(crossprod(Z2) - diag(r2))) > res){
    message = "The quantified constraints are provided by the orthogonality of Z2."
    res = max(abs(crossprod(Z2) - diag(r2)))
  }
  
  # Check whether the product of U1 and U2 is diagonal
  if (r0 > 1){
    U12 = crossprod(U1,U2)
    temp_res = -Inf
    for (i in 1:r0){
      temp = max(abs(U12[i,-i]))
      if (temp > temp_res){
        temp_res = temp
      }
    }
    if (temp_res > res){
      message = "The quantified constraints are provided by the non-diagonal entry of crossprod(U1,U2)."
      res = temp_res
    }
  }
  return(list(res = res, msg = message))
}

#' Function that calculates principal angles between column spaces of two matrices
#'
#' @param X The first matrix.
#' @param Y The second matrix.
#' @param tol Tolerence, default is the square root of machine precision.
#'
#' @return A list that contains the following:
#' \item{angle}{A vector of principal angles with increasing order.}
#' \item{cos_angle}{A vector of cosine principal angles}
#' \item{principal_vector1}{Principal vectors of matrix \code{X}}
#' \item{principal_vector2}{Principal vectors of matrix \code{Y}}
#'
#' @examples
#' X = matrix(c(1,1,1,1,1,0),nrow = 3, ncol = 2)
#' Y = matrix(c(1,1,1,2,1,0),nrow = 3, ncol = 2)
#' angle_cal(X,Y)
#' 
angle_cal <- function(X, Y, tol = .Machine$double.eps^0.5){
  X = as.matrix(X)
  Y = as.matrix(Y)
  Xnorm = svd(X)$u
  Ynorm = svd(Y)$u
  M = crossprod(Xnorm,Ynorm)
  # Extreme case when both X and Y only contain one number
  if (dim(M)[1] == 1 && dim(M)[2] == 1){
    cos_angle = abs(M)
    principal_angle = NA
    if (cos_angle >= 1){principal_angle = 0}
    if (cos_angle <= 0){principal_angle = pi/2}
    if (cos_angle > 0 && cos_angle < 1){
      principal_angle = acos(cos_angle)
    }
    principal_mat1 = Xnorm 
    principal_mat2 = Ynorm 
    return(list("angle" = principal_angle, "cos_angle" = cos_angle,
                "principal_vector1" = principal_mat1, 
                "principal_vector2" = principal_mat2))
  }
  # Normal case when X and Y are matrices (data frames)
  else{
    svd_result = svd(M)
    cos_angle = svd_result$d
    l = length(cos_angle)
    principal_angle = rep(NA, l)
    for (i in 1:l){
      if (cos_angle[i] >= 1){principal_angle[i] = 0}
      if (cos_angle[i] <= 0){principal_angle[i] = pi/2}
      if (cos_angle[i] > 0 && cos_angle[i] < 1){
        principal_angle[i] = acos(cos_angle[i])
      }
    }
    principal_mat1 = Xnorm %*% svd_result$u
    principal_mat2 = Ynorm %*% svd_result$v
    return(list("angle" = principal_angle, "cos_angle" = cos_angle,
                "principal_vector1" = principal_mat1, 
                "principal_vector2" = principal_mat2))
  }
}

# Calculate the projection matrix of X.
projection <- function(X){
  return(X %*% pinv(X))
}