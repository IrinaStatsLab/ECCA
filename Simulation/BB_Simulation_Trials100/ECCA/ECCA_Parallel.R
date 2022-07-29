# The ECCA result for this setting is run on the server with parallel computing. Replace the path to be your own path.

rm(list = ls())
source('../../../exp_family.R')
source('../../../exponentialCCA_v2.R')
source('../../../optU.R')
source('../../../SOC.R')
source('../../../newton_method_new.R')
source('../../../AnalyticalUpdate.R')

# Load the generated data
load('../Data.RData')
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

library(foreach)
library(doParallel)
cl = makeCluster(28)
registerDoParallel(cl)

ones = rep(1, n)

output <- foreach (i = 1:nrep) %dopar% {
  library(pracma)
  X1 = X1_true[[i]]
  X2 = X2_true[[i]]
  X1 = modifyBinomial(X1,trials)
  X2 = modifyBinomial(X2,trials)
  
  temp_res = ecca_v2(X1, X2, r0, r1, r2, family1 = family1, family2 = family2, gamma = 1000, trials = trials)
  theta1 = outer(ones, temp_res$mu1) + tcrossprod(temp_res$U1, temp_res$V1) + tcrossprod(temp_res$Z1, temp_res$A1)
  theta2 = outer(ones, temp_res$mu2) + tcrossprod(temp_res$U2, temp_res$V2) + tcrossprod(temp_res$Z2, temp_res$A2)

  list(theta1_est = theta1,
       theta2_est = theta2,
       U1_est = temp_res$U1,
       U2_est = temp_res$U2,
       Z1_est = temp_res$Z1,
       Z2_est = temp_res$Z2,
       mu1_est = temp_res$mu1,
       mu2_est = temp_res$mu2,
       V1_est = temp_res$V1,
       V2_est = temp_res$V2,
       A1_est = temp_res$A1,
       A2_est = temp_res$A2,
       loss = temp_res$loss_trace,
       constraints = temp_res$constraints)
}

stopCluster(cl)
save(output, file = "ECCA_result.RData")