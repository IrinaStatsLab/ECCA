rm(list = ls())
source("Simulation/data_generator.R")
source("Function/exp_family.R")
set.seed(37)
# Set parameters
n = 10
p1 = 8
p2 = 7
p = c(p1, p2)
family1 = 'gaussian'
family2 = 'binomial'
lambdas = c(1, 0.7) # covariance between joint structures
trials = 100

# Set ranks (joint and individuals)
r0 = 2
r1 = 1
r2 = 2
SNR = 5

# Generate data
sampledata = dataGenerator(r0 = r0, r1 = r1, r2 = r2, n = n, p = p, lambda = lambdas, family1, family2, SNR = SNR, size = trials, main_effect = TRUE)
# Get true signal
theta1_true = sampledata$theta1
theta2_true = sampledata$theta2
Z1_true = sampledata$Z1
Z2_true = sampledata$Z2
U1_true = sampledata$U1
U2_true = sampledata$U2

# Sanity check
print(crossprod(U1_true, U2_true))

# Source all necessary functions to apply ECCA
source('Function/exponentialCCA.R')
source('Function/optU.R')
source('Function/SOC.R')
source('Function/newton_method.R')
source('Function/AnalyticalUpdate.R')
# Apply ECCA
ECCA_result = ecca_v2(x1 = sampledata$X1, x2 = sampledata$X2, r0 = r0, r1 = r1, r2 = r2, family1 = family1, family2 = family2, trials = trials)

# Record ECCA results
U1_ecca = ECCA_result$U1
U2_ecca = ECCA_result$U2
Z1_ecca = ECCA_result$Z1
Z2_ecca = ECCA_result$Z2 

one_vec = rep(1, n)

# Get estimated natural parameters
theta1_ecca = tcrossprod(one_vec, ECCA_result$mu1) + tcrossprod(U1_ecca, ECCA_result$V1) + tcrossprod(Z1_ecca, ECCA_result$A1)
theta2_ecca = tcrossprod(one_vec, ECCA_result$mu2) + tcrossprod(U2_ecca, ECCA_result$V2) + tcrossprod(Z2_ecca, ECCA_result$A2)

# Check the relative error between estimated signal of DMMD with true signal 
Fnorm(theta1_ecca - theta1_true)^2/Fnorm(theta1_true)^2
Fnorm(theta2_ecca - theta2_true)^2/Fnorm(theta2_true)^2

# Check the chordal distance between estimated joint space
1/sqrt(2)*Fnorm(projection(U1_ecca) - projection(U1_true))
1/sqrt(2)*Fnorm(projection(U2_ecca) - projection(U2_true))

# Check the estimated canonical correlations
angles = angle_cal(U1_ecca, U2_ecca)$angle
cos(angles)
