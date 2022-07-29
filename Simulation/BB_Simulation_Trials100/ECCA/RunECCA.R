# This file is only for illustration purpose. The actual result is run on the server with parallel computing. See 'ECCA_Parallel.R' file.
rm(list = ls())
source('Function/exp_family.R')
source('Function/exponentialCCA.R')
source('Function/optU.R')
source('Function/SOC.R')
source('Function/newton_method.R')
source('Function/AnalyticalUpdate.R')

# Load the generated data
load('Simulation/BB_Simulation_Trials100/Data.RData')
# Set parameters
n = 50
p1 = 30
p2 = 20
p = c(p1, p2)
family1 = 'binomial'
family2 = 'binomial'

# Set ranks (joint and individuals)
r0 = 3
r1 = 4
r2 = 3
trials = 100

nrep = 100
set.seed(37)

# Initialize the saving lists.
theta1_est = list()
theta2_est = list()
U1_est = list()
U2_est = list()
Z1_est = list()
Z2_est = list()
mu1_est = list()
mu2_est = list()
V1_est = list()
V2_est = list()
A1_est = list()
A2_est = list()
loss = list()
constraints = list()

ones = rep(1, n)
# Take some time to run.
for (i in 1:nrep){
  print(i)
  X1 = X1_true[[i]]
  X2 = X2_true[[i]]
  X1 = modifyBinomial(X1,trials)
  X2 = modifyBinomial(X2,trials)
  
  temp_res = ecca_v2(X1, X2, r0, r1, r2, family1 = family1, family2 = family2, gamma = 1000, trials = trials)
  theta1 = outer(ones, temp_res$mu1) + tcrossprod(temp_res$U1, temp_res$V1) + tcrossprod(temp_res$Z1, temp_res$A1)
  theta2 = outer(ones, temp_res$mu2) + tcrossprod(temp_res$U2, temp_res$V2) + tcrossprod(temp_res$Z2, temp_res$A2)
  theta1_est = append(theta1_est, list(theta1))
  theta2_est = append(theta2_est, list(theta2))

  U1_est = append(U1_est, list(temp_res$U1))
  U2_est = append(U2_est, list(temp_res$U2))
  Z1_est = append(Z1_est, list(temp_res$Z1))
  Z2_est = append(Z2_est, list(temp_res$Z2))
  mu1_est = append(mu1_est, list(temp_res$mu1))
  mu2_est = append(mu2_est, list(temp_res$mu2))
  V1_est = append(V1_est, list(temp_res$V1))
  V2_est = append(V2_est, list(temp_res$V2))
  A1_est = append(A1_est, list(temp_res$A1))
  A2_est = append(A2_est, list(temp_res$A2))
  loss = append(loss, list(temp_res$loss_trace))
  constraints = append(constraints, list(temp_res$constraints))
}

save(theta1_est, theta2_est, U1_est, U2_est, Z1_est, Z2_est, mu1_est, mu2_est, 
     V1_est, V2_est, A1_est, A2_est, loss, constraints, file = "Simulation/BB_Simulation_Trials100/ECCA_result.RData")