rm(list = ls())
# Run on the cluster
source('Simulation/data_generator.R')
source('Function/exp_family.R')

# Set parameters
n = 50
p1 = 30
p2 = 20
p = c(p1, p2)
family1 = 'gaussian'
family2 = 'binomial'
lambdas = c(1, 0.9, 0.7) # covariance between joint structures
trials = 100

# Set ranks (joint and individuals)
r0 = 3
r1 = 4
r2 = 3
SNR = 5

nrep = 100
set.seed(37)

# Initialize the saving lists.

X1_true = list()
X2_true = list()
theta1_true = list()
theta2_true = list()
U1_true = list()
U2_true = list()
Z1_true = list()
Z2_true = list()
mu1_true = list()
mu2_true = list()
V1_true = list()
V2_true = list()
A1_true = list()
A2_true = list()

for (i in 1:nrep){
  sampledata = dataGenerator(r0, r1, r2, n = n, p = p, lambda = lambdas, family1, family2, SNR = SNR, size = trials, main_effect = TRUE)
  
  X1_true = append(X1_true, list(sampledata$X1))
  X2_true = append(X2_true, list(sampledata$X2))
  theta1_true = append(theta1_true, list(sampledata$theta1))
  theta2_true = append(theta2_true, list(sampledata$theta2))
  U1_true = append(U1_true, list(sampledata$U1))
  U2_true = append(U2_true, list(sampledata$U2))
  Z1_true = append(Z1_true, list(sampledata$Z1))
  Z2_true = append(Z2_true, list(sampledata$Z2))
  mu1_true = append(mu1_true, list(sampledata$mu1))
  mu2_true = append(mu2_true, list(sampledata$mu2))
  V1_true = append(V1_true, list(sampledata$V1))
  V2_true = append(V2_true, list(sampledata$V2))
  A1_true = append(A1_true, list(sampledata$A1))
  A2_true = append(A2_true, list(sampledata$A2))
}

save(X1_true, X2_true, theta1_true, theta2_true, U1_true, U2_true, Z1_true, Z2_true, 
     mu1_true, mu2_true, V1_true, V2_true, A1_true, A2_true, file = "Simulation/GB_Simulation_SNR5/Data.RData")