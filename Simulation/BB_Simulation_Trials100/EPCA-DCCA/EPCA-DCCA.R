rm(list = ls())
library(generalizedPCA)
# Run python code
library("reticulate")

# Check the version of Python.
py_config()

# Source the dcca function
source_python("Simulation/dcca.py")
source('Function/exp_family.R')
# Load the generated data
load('Simulation/BB_Simulation_Trials100/Data.RData')

set.seed(37)
# Set parameters
n = 50
p1 = 30
p2 = 20
p = c(p1, p2)
family1 = 'binomial'
family2 = 'binomial'
trials = 100

# Set ranks (joint and individuals)
r0 = 3
r1 = 4
r2 = 3

r0_rank1 = 1
r1_rank1 = 6
r2_rank1 = 5

ones = rep(1, n)
nrep = 100
set.seed(37)

# Initialize the saving lists.
res_DCCA_list = list()
res_EPCA_DCCA_list = list()

res_DCCA_list_rank1 = list()
res_EPCA_DCCA_list_rank1 = list()

# Save the main effect for DCCA:
res_DCCA_mu1_list = list()
res_DCCA_mu2_list = list()

# Save the main effect for EPCA-DCCA:
res_EPCA_DCCA_mu1_list = list()
res_EPCA_DCCA_mu2_list = list()


for (i in 1:nrep){
  X1 = X1_true[[i]]
  X2 = X2_true[[i]]
  X1 = modifyBinomial(X1, trials)
  X2 = modifyBinomial(X2, trials)
  
  para1 = generalizedPCA(X1, k = r0 + r1, family = family1)
  para2 = generalizedPCA(X2, k = r0 + r2, family = family2)
  
  theta1_est = saturated_para(X1, family1, trials = trials)
  theta2_est = saturated_para(X2, family2, trials = trials)
  
  # update main effect:
  res_EPCA_DCCA_mu1_list = append(res_EPCA_DCCA_mu1_list, list(para1$mu))
  res_EPCA_DCCA_mu2_list = append(res_EPCA_DCCA_mu2_list, list(para2$mu))
  
  res_DCCA_mu1_list = append(res_DCCA_mu1_list, list(colMeans(theta1_est)))
  res_DCCA_mu2_list = append(res_DCCA_mu2_list, list(colMeans(theta2_est)))
  
  theta1_epca = (theta1_est - tcrossprod(ones, para1$mu)) %*% tcrossprod(para1$U)
  theta2_epca = (theta2_est - tcrossprod(ones, para2$mu)) %*% tcrossprod(para2$U)
  
  theta1 = theta1_est - tcrossprod(ones, colMeans(theta1_est))
  theta2 = theta2_est - tcrossprod(ones, colMeans(theta2_est))
  
  res_DCCA = dCCA(t(theta1), t(theta2), r_1 = r0 + r1, r_2 = r0 + r2, r_12 = r0)
  res_EPCA_DCCA = dCCA(t(theta1_epca), t(theta2_epca), r_1 = r0 + r1, r_2 = r0 + r2, r_12 = r0)
  
  res_DCCA_list = append(res_DCCA_list, list(res_DCCA))
  res_EPCA_DCCA_list = append(res_EPCA_DCCA_list, list(res_EPCA_DCCA))
  
  res_DCCA_rank1 = dCCA(t(theta1), t(theta2), r_1 = r0 + r1, r_2 = r0 + r2, r_12 = 1)
  res_EPCA_DCCA_rank1 = dCCA(t(theta1_epca), t(theta2_epca), r_1 = r0 + r1, r_2 = r0 + r2, r_12 = 1)
  
  res_DCCA_list_rank1 = append(res_DCCA_list_rank1, list(res_DCCA_rank1))
  res_EPCA_DCCA_list_rank1 = append(res_EPCA_DCCA_list_rank1, list(res_EPCA_DCCA_rank1))
}
save(res_DCCA_list, res_EPCA_DCCA_list, res_DCCA_list_rank1, res_EPCA_DCCA_list_rank1, 
     res_DCCA_mu1_list, res_DCCA_mu2_list, res_EPCA_DCCA_mu1_list, res_EPCA_DCCA_mu2_list, file = "Simulation/BB_Simulation_Trials/DCCA_result.RData")