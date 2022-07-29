#----------------------------------
# Dongbang Yuan
# Contact: yuandb09@gmail.com
#------------------------------

# This is the analytical update for score matrices Z and U to replace Newton's method in SOC update. 
# Notice the optimization problem is Augmented Lagrangian.
# This can be used when either exponential family is Gaussian.
# Refer to Appendix B.2.3 for more details. The notation follows.
# For example, X1 is Gaussian, we can use this function to update Z1, although we still need to follow SOC update to get P1, B1.
Analyticalupdate <- function(X, mu = NULL, U, V, B, P, A, gamma, main_effect = TRUE){
  n = dim(X)[1]
  r = dim(A)[2]
  if (main_effect){
    if (is.null(mu)){
      stop('Please provide the value for the mean vector')
    }
    one_vector = rep(1, n)
    result = ((X - tcrossprod(one_vector, mu) - tcrossprod(U, V)) %*% A - gamma * (B - P)) %*% solve(crossprod(A) + gamma * diag(r))
  }
  else{
    result = ((X - tcrossprod(U, V)) %*% A - gamma * (B - P)) %*% solve(crossprod(A) + gamma * diag(r))
  }
  return(result)
}

# This is the analytical update for updating loading matrices.
# This is used only when the exponential family is Gaussian.
# Refer to Appendix B.1.1 for more details. The notation follows.
# For example, X1 is Gaussian, we can use this function to update mu1, V1 and A1.
require(pracma)
MainAnalyticalupdate <- function(X, S, r_joint, r_ind, main_effect = TRUE){
  n = dim(X)[1]
  if (dim(X)[1] != dim(S)[1]){
    stop("The dimension of X does not agree with the dimension of S.")
  }
  whole_loading = t(pinv(S) %*% X)
  if (main_effect){
    mu = whole_loading[,1]
    if (r_joint == 0){
      V = matrix(0, nrow = n, ncol = 1)
    }
    else{
      V = whole_loading[,2:(1+r_joint),drop = F]
    }
    if (r_ind == 0){
      A = matrix(0, nrow = n, ncol = 1)
    }
    else{
      A = whole_loading[,(2+r_joint):(1+r_joint+r_ind),drop = F]
    }
    return(list("mu" = mu, "V" = V, "A" = A))
  }
  else{
    if (r_joint == 0){
      V = matrix(0, nrow = n, ncol = 1)
    }
    else{
      V = whole_loading[,1:r_joint,drop = F]
    }
    if (r_ind == 0){
      A = matrix(0, nrow = n, ncol = 1)
    }
    else{
      A = whole_loading[,(1+r_joint):(r_joint+r_ind),drop = F]
    }
    return(list("V" = V, "A" = A))
  }
}

# This is the analytical update for Z to replace SOC completely.
# This works when the families are BOTH Gaussian, i.e., X1 and X2 are both Gaussian.
# Refer to Appendix B.2.1 for more details. The notation follows.
AnalyticalSOC_Z <- function(x1, x2, mu1 = NULL, mu2 = NULL, U1, U2, V1, V2, A1, A2, main_effect = TRUE){
  n = dim(x1)[1]
  p1 = dim(x1)[2]
  p2 = dim(x2)[2]
  r1 = dim(A1)[2]
  r2 = dim(A2)[2]
  if (main_effect){
    one_vector = rep(1, n)
    U = cbind(one_vector,U1,U2)
  }
  else{
    U = cbind(U1,U2)
  }
  projU_ortho = diag(n) - U %*% pinv(U)
  zero_matrix1 = matrix(rep(0,p1*r2), nrow = p1)
  zero_matrix2 = matrix(rep(0,p2*r1), nrow = p2)
  A = rbind(cbind(A1, zero_matrix1), cbind(zero_matrix2, A2))
  if (main_effect){
    if (is.null(mu1)|is.null(mu2)){
      stop('Please provide the value for the mean vectors')
    }
    Y1 = x1 - tcrossprod(one_vector, mu1) - tcrossprod(U1, V1)
    Y2 = x2 - tcrossprod(one_vector, mu2) - tcrossprod(U2, V2)
  }
  else{
    Y1 = x1 - tcrossprod(U1, V1)
    Y2 = x2 - tcrossprod(U2, V2)
  }
  svd_result = svd(projU_ortho %*% cbind(Y1, Y2) %*% A)
  result = tcrossprod(svd_result$u, svd_result$v)
  if (r1 == 0){
    Z1 = matrix(rep(0, n), nrow = n)
  }
  else{
    Z1 = result[,1:r1, drop = F]
  }
  if (r2 == 0){
    Z2 = matrix(rep(0, n), nrow = n)
  }
  else{
    Z2 = result[,(r1+1):(r1+r2), drop = F]
  }
  return(list(Z1 = Z1, Z2 = Z2))
}

# This is the analytical update for U to replace SOC completely.
# It works when either family is Gaussian.
# When first is TRUE, it means we update U1; if False, we update U2.
# For example, X1 is Gaussian, we can use this function to update U1. (first = TRUE).
# Refer to Appendix B.3.1 for more details. The notation follows.
AnalyticalSOC_U <- function(x, mu = NULL, Z1, Z2, V, A, first = TRUE, main_effect = TRUE){
  n = dim(x)[1]
  if (main_effect){
    if (is.null(mu)){
      stop('Please provide the value of the mean vector.')
    }
    one_vector = rep(1, n)
    Z = cbind(one_vector, Z1, Z2)
    projZ_ortho = diag(n) - Z %*% pinv(Z)
    if (first){
      B = x - tcrossprod(one_vector,mu) - tcrossprod(Z1,A)
    }
    else{
      B = x - tcrossprod(one_vector,mu) - tcrossprod(Z2,A)
    }
    svd_result = svd(projZ_ortho %*% B %*% V)
    U = tcrossprod(svd_result$u, svd_result$v)
  }
  else{
    Z = cbind(Z1,Z2)
    # Here, to get the projection matrix onto the orthogonal complement of Z, we could use I - ZZ^T,
    projZ_ortho = diag(n) - Z %*% pinv(Z)
    if (first){
      B = x - tcrossprod(Z1,A)
    }
    else{
      B = x - tcrossprod(Z2,A)
    }
    svd_result = svd(projZ_ortho %*% B %*% V)
    U = tcrossprod(svd_result$u, svd_result$v)
  }
  return(U)
}