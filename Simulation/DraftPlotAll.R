rm(list = ls())
source('Function/exp_family.R')
n = 50
p1 = 30
p2 = 20
p = c(p1, p2)
lambdas = c(1, 0.9, 0.7) # covariance between joint structures
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

############################## Gaussian - Gaussian ##############################
source('Simulation/GG_Simulation_SNR5/GAS_GG_results.RData')
load('Simulation/GG_Simulation_SNR5/DCCA_result.RData')
load('Simulation/GG_Simulation_SNR5/ECCA_result.RData')
load('Simulation/GG_Simulation_SNR5/Data.RData')

# Initialize the signal error vectors.
ECCA_signal_error1_gg = rep(NA, nrep)
ECCA_signal_error2_gg = rep(NA, nrep)

DCCA_signal_error1_gg = rep(NA, nrep)
DCCA_signal_error2_gg = rep(NA, nrep)

EPCA_DCCA_signal_error1_gg = rep(NA, nrep)
EPCA_DCCA_signal_error2_gg = rep(NA, nrep)

GAS_signal_error1_gg = rep(NA, nrep)
GAS_signal_error2_gg = rep(NA, nrep)

DCCA_signal_error1_rank1_gg = rep(NA, nrep)
DCCA_signal_error2_rank1_gg = rep(NA, nrep)

EPCA_DCCA_signal_error1_rank1_gg = rep(NA, nrep)
EPCA_DCCA_signal_error2_rank1_gg = rep(NA, nrep)

GAS_signal_error1_rank1_gg = rep(NA, nrep)
GAS_signal_error2_rank1_gg = rep(NA, nrep)

###################### Compare the correlated signals ######################
# The idea is to get the first three canonical variables from estimated centered natural parameter matrix from GAS, DCCA, EPCA-DCCA 
# and compare it with the true 'correlated' signals.

# Initialize joint error vector.
ECCA_joint_error1_gg = rep(NA, nrep)
ECCA_joint_error2_gg = rep(NA, nrep)

DCCA_joint_error1_gg = rep(NA, nrep)
DCCA_joint_error2_gg = rep(NA, nrep)

EPCA_DCCA_joint_error1_gg = rep(NA, nrep)
EPCA_DCCA_joint_error2_gg = rep(NA, nrep)

GAS_joint_error1_gg = rep(NA, nrep)
GAS_joint_error2_gg = rep(NA, nrep)

# For joint rank 1 result:
DCCA_joint_error1_rank1_gg = rep(NA, nrep)
DCCA_joint_error2_rank1_gg = rep(NA, nrep)

EPCA_DCCA_joint_error1_rank1_gg = rep(NA, nrep)
EPCA_DCCA_joint_error2_rank1_gg = rep(NA, nrep)

GAS_joint_error1_rank1_gg = rep(NA, nrep)
GAS_joint_error2_rank1_gg = rep(NA, nrep)


###################### ECCA ##########################
for (i in 1:nrep){
  ECCA_signal_error1_gg[i] = Fnorm(theta1_est[[i]] - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  ECCA_signal_error2_gg[i] = Fnorm(theta2_est[[i]] - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  joint1_ecca = U1_est[[i]]
  joint2_ecca = U2_est[[i]]
  
  ECCA_joint_error1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_ecca) - projection(joint1_true))
  ECCA_joint_error2_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_ecca) - projection(joint2_true))
}

###################### GAS ##########################
for (i in 1:nrep){
  U0 = U0_GAS[,((i-1)*r0 + 1):(i*r0)]
  U1 = U1_GAS[,((i-1)*r1 + 1):(i*r1)]
  U2 = U2_GAS[,((i-1)*r2 + 1):(i*r2)]
  
  V1 = V0_GAS[1:p1,((i-1)*r0 + 1):(i*r0)]
  V2 = V0_GAS[(p1+1):(p1+p2),((i-1)*r0 + 1):(i*r0)]
  
  A1 = A1_GAS[,((i-1)*r1 + 1):(i*r1)]
  A2 = A2_GAS[,((i-1)*r2 + 1):(i*r2)]
  
  mu1 = Mu0_GAS[1:p1,i]
  mu2 = Mu0_GAS[(p1+1):(p1+p2),i]
  
  theta1_gas = tcrossprod(ones, mu1) + tcrossprod(U0, V1) + tcrossprod(U1, A1)
  theta2_gas = tcrossprod(ones, mu2) + tcrossprod(U0, V2) + tcrossprod(U2, A2)
  
  GAS_signal_error1_gg[i] = Fnorm(theta1_gas - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_gg[i] = Fnorm(theta2_gas - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  U0_rank1 = U0_GAS_rank1[,((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  U1_rank1 = U1_GAS_rank1[,((i-1)*r1_rank1 + 1):(i*r1_rank1)]
  U2_rank1 = U2_GAS_rank1[,((i-1)*r2_rank1 + 1):(i*r2_rank1)]
  
  V1_rank1 = V0_GAS_rank1[1:p1,((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  V2_rank1 = V0_GAS_rank1[(p1+1):(p1+p2),((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  
  A1_rank1 = A1_GAS_rank1[,((i-1)*r1_rank1 + 1):(i*r1_rank1)]
  A2_rank1 = A2_GAS_rank1[,((i-1)*r2_rank1 + 1):(i*r2_rank1)]
  
  mu1_rank1 = Mu0_GAS_rank1[1:p1,i]
  mu2_rank1 = Mu0_GAS_rank1[(p1+1):(p1+p2),i]
  
  theta1_gas_rank1 = tcrossprod(ones, mu1_rank1) + tcrossprod(U0_rank1, V1_rank1) + tcrossprod(U1_rank1, A1_rank1)
  theta2_gas_rank1 = tcrossprod(ones, mu2_rank1) + tcrossprod(U0_rank1, V2_rank1) + tcrossprod(U2_rank1, A2_rank1)
  
  GAS_signal_error1_rank1_gg[i] = Fnorm(theta1_gas_rank1 - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_rank1_gg[i] = Fnorm(theta2_gas_rank1 - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  theta1_gas_center = tcrossprod(U0, V1) + tcrossprod(U1, A1)
  theta2_gas_center = tcrossprod(U0, V2) + tcrossprod(U2, A2)
  
  theta1_gas_center_rank1 = tcrossprod(U0_rank1, V1_rank1) + tcrossprod(U1_rank1, A1_rank1)
  theta2_gas_center_rank1 = tcrossprod(U0_rank1, V2_rank1) + tcrossprod(U2_rank1, A2_rank1)
  
  col1 = svd(theta1_gas_center)$u[,1:(r0+r1)]
  col2 = svd(theta2_gas_center)$u[,1:(r0+r2)]
  angle_result = angle_cal(col1, col2)
  
  col1_rank1 = svd(theta1_gas_center_rank1)$u[,1:(r0+r1)]
  col2_rank1 = svd(theta2_gas_center_rank1)$u[,1:(r0+r2)]
  angle_result_rank1 = angle_cal(col1_rank1, col2_rank1)
  
  # Get GAS canonical variables
  joint1_GAS = angle_result$principal_vector1[,1:r0]
  joint2_GAS = angle_result$principal_vector2[,1:r0]
  
  joint1_GAS_rank1 = angle_result_rank1$principal_vector1[,1:r0]
  joint2_GAS_rank1 = angle_result_rank1$principal_vector2[,1:r0]
  
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  GAS_joint_error1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS) - projection(joint1_true))
  GAS_joint_error2_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS) - projection(joint2_true))
  
  GAS_joint_error1_rank1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS_rank1) - projection(joint1_true))
  GAS_joint_error2_rank1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS_rank1) - projection(joint2_true))
}

###################### DCCA ##########################
for (i in 1:nrep){
  DCCA_signal_error1_gg[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_gg[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_gg[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_gg[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  DCCA_signal_error1_rank1_gg[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_rank1_gg[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_rank1_gg[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_rank1_gg[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
}

for (i in 1:nrep){
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  theta1_dcca_center = t(res_DCCA_list[[i]][[1]])
  theta2_dcca_center = t(res_DCCA_list[[i]][[2]])
  
  theta1_epca_dcca_center = t(res_EPCA_DCCA_list[[i]][[1]])
  theta2_epca_dcca_center = t(res_EPCA_DCCA_list[[i]][[2]])
  
  theta1_dcca_center_rank1 = t(res_DCCA_list_rank1[[i]][[1]])
  theta2_dcca_center_rank1 = t(res_DCCA_list_rank1[[i]][[2]])
  
  theta1_epca_dcca_center_rank1 = t(res_EPCA_DCCA_list_rank1[[i]][[1]])
  theta2_epca_dcca_center_rank1 = t(res_EPCA_DCCA_list_rank1[[i]][[2]])
  
  # DCCA joint
  col1_dcca = svd(theta1_dcca_center)$u[,1:(r0+r1)]
  col2_dcca = svd(theta2_dcca_center)$u[,1:(r0+r2)]
  angle_result_dcca = angle_cal(col1_dcca, col2_dcca)
  
  col1_dcca_rank1 = svd(theta1_dcca_center_rank1)$u[,1:(r0+r1)]
  col2_dcca_rank1 = svd(theta2_dcca_center_rank1)$u[,1:(r0+r2)]
  angle_result_dcca_rank1 = angle_cal(col1_dcca_rank1, col2_dcca_rank1)
  
  # EPCA-DCCA joint 
  col1_epca_dcca = svd(theta1_epca_dcca_center)$u[,1:(r0+r1)]
  col2_epca_dcca = svd(theta2_epca_dcca_center)$u[,1:(r0+r2)]
  angle_result_epca_dcca = angle_cal(col1_epca_dcca, col2_epca_dcca)
  
  col1_epca_dcca_rank1 = svd(theta1_epca_dcca_center_rank1)$u[,1:(r0+r1)]
  col2_epca_dcca_rank1 = svd(theta2_epca_dcca_center_rank1)$u[,1:(r0+r2)]
  angle_result_epca_dcca_rank1 = angle_cal(col1_epca_dcca_rank1, col2_epca_dcca_rank1)
  
  # Get DCCA canonical variables
  joint1_dcca = angle_result_dcca$principal_vector1[,1:r0]
  joint2_dcca = angle_result_dcca$principal_vector2[,1:r0]
  
  joint1_dcca_rank1 = angle_result_dcca_rank1$principal_vector1[,1:r0]
  joint2_dcca_rank1 = angle_result_dcca_rank1$principal_vector2[,1:r0]
  
  # Get EPCA-DCCA canonical variables
  joint1_epca_dcca = angle_result_epca_dcca$principal_vector1[,1:r0]
  joint2_epca_dcca = angle_result_epca_dcca$principal_vector2[,1:r0]
  
  joint1_epca_dcca_rank1 = angle_result_epca_dcca_rank1$principal_vector1[,1:r0]
  joint2_epca_dcca_rank1 = angle_result_epca_dcca_rank1$principal_vector2[,1:r0]
  
  #### Calculate the chordal distance between subspaces.
  DCCA_joint_error1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca) - projection(joint1_true))
  DCCA_joint_error2_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca) - projection(joint2_true))
  
  DCCA_joint_error1_rank1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca_rank1) - projection(joint1_true))
  DCCA_joint_error2_rank1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca_rank1) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca) - projection(joint1_true))
  EPCA_DCCA_joint_error2_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_rank1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca_rank1) - projection(joint1_true))
  EPCA_DCCA_joint_error2_rank1_gg[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca_rank1) - projection(joint2_true))
}

############################## Gaussian - Binomial ##############################
source('Simulation/GB_Simulation_SNR5/GAS_GB_results.RData')
load('Simulation/GB_Simulation_SNR5/DCCA_result.RData')
load('Simulation/GB_Simulation_SNR5/ECCA_result.RData')
load('Simulation/GB_Simulation_SNR5/Data.RData')

trials = 100

# Initialize the signal error vectors.
ECCA_signal_error1_gb = rep(NA, nrep)
ECCA_signal_error2_gb = rep(NA, nrep)

DCCA_signal_error1_gb = rep(NA, nrep)
DCCA_signal_error2_gb = rep(NA, nrep)

EPCA_DCCA_signal_error1_gb = rep(NA, nrep)
EPCA_DCCA_signal_error2_gb = rep(NA, nrep)

GAS_signal_error1_gb = rep(NA, nrep)
GAS_signal_error2_gb = rep(NA, nrep)

DCCA_signal_error1_rank1_gb = rep(NA, nrep)
DCCA_signal_error2_rank1_gb = rep(NA, nrep)

EPCA_DCCA_signal_error1_rank1_gb = rep(NA, nrep)
EPCA_DCCA_signal_error2_rank1_gb = rep(NA, nrep)

GAS_signal_error1_rank1_gb = rep(NA, nrep)
GAS_signal_error2_rank1_gb = rep(NA, nrep)

###################### Compare the correlated signals ######################
# The idea is to get the first three canonical variables from estimated centered natural parameter matrix from GAS, DCCA, EPCA-DCCA 
# and compare it with the true 'correlated' signals.

# Initialize joint error vector.
ECCA_joint_error1_gb = rep(NA, nrep)
ECCA_joint_error2_gb = rep(NA, nrep)

DCCA_joint_error1_gb = rep(NA, nrep)
DCCA_joint_error2_gb = rep(NA, nrep)

EPCA_DCCA_joint_error1_gb = rep(NA, nrep)
EPCA_DCCA_joint_error2_gb = rep(NA, nrep)

GAS_joint_error1_gb = rep(NA, nrep)
GAS_joint_error2_gb = rep(NA, nrep)

# For joint rank 1 result:
DCCA_joint_error1_rank1_gb = rep(NA, nrep)
DCCA_joint_error2_rank1_gb = rep(NA, nrep)

EPCA_DCCA_joint_error1_rank1_gb = rep(NA, nrep)
EPCA_DCCA_joint_error2_rank1_gb = rep(NA, nrep)

GAS_joint_error1_rank1_gb = rep(NA, nrep)
GAS_joint_error2_rank1_gb = rep(NA, nrep)

###################### ECCA ##########################
for (i in 1:nrep){
  ECCA_signal_error1_gb[i] = Fnorm(output[[i]]$theta1_est - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  ECCA_signal_error2_gb[i] = Fnorm(output[[i]]$theta2_est - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  joint1_ecca = output[[i]]$U1_est
  joint2_ecca = output[[i]]$U2_est
  
  ECCA_joint_error1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_ecca) - projection(joint1_true))
  ECCA_joint_error2_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_ecca) - projection(joint2_true))
}

###################### GAS ##########################
for (i in 1:nrep){
  U0 = U0_GAS[,((i-1)*r0 + 1):(i*r0)]
  U1 = U1_GAS[,((i-1)*r1 + 1):(i*r1)]
  U2 = U2_GAS[,((i-1)*r2 + 1):(i*r2)]
  
  V1 = V0_GAS[1:p1,((i-1)*r0 + 1):(i*r0)]
  V2 = V0_GAS[(p1+1):(p1+p2),((i-1)*r0 + 1):(i*r0)]
  
  A1 = A1_GAS[,((i-1)*r1 + 1):(i*r1)]
  A2 = A2_GAS[,((i-1)*r2 + 1):(i*r2)]
  
  mu1 = Mu0_GAS[1:p1,i]
  mu2 = Mu0_GAS[(p1+1):(p1+p2),i]
  
  theta1_gas = tcrossprod(ones, mu1) + tcrossprod(U0, V1) + tcrossprod(U1, A1)
  theta2_gas = trials * (tcrossprod(ones, mu2) + tcrossprod(U0, V2) + tcrossprod(U2, A2))
  
  GAS_signal_error1_gb[i] = Fnorm(theta1_gas - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_gb[i] = Fnorm(theta2_gas - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  U0_rank1 = U0_GAS_rank1[,((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  U1_rank1 = U1_GAS_rank1[,((i-1)*r1_rank1 + 1):(i*r1_rank1)]
  U2_rank1 = U2_GAS_rank1[,((i-1)*r2_rank1 + 1):(i*r2_rank1)]
  
  V1_rank1 = V0_GAS_rank1[1:p1,((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  V2_rank1 = V0_GAS_rank1[(p1+1):(p1+p2),((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  
  A1_rank1 = A1_GAS_rank1[,((i-1)*r1_rank1 + 1):(i*r1_rank1)]
  A2_rank1 = A2_GAS_rank1[,((i-1)*r2_rank1 + 1):(i*r2_rank1)]
  
  mu1_rank1 = Mu0_GAS_rank1[1:p1,i]
  mu2_rank1 = Mu0_GAS_rank1[(p1+1):(p1+p2),i]
  
  theta1_gas_rank1 = tcrossprod(ones, mu1_rank1) + tcrossprod(U0_rank1, V1_rank1) + tcrossprod(U1_rank1, A1_rank1)
  theta2_gas_rank1 = trials * (tcrossprod(ones, mu2_rank1) + tcrossprod(U0_rank1, V2_rank1) + tcrossprod(U2_rank1, A2_rank1))
  
  GAS_signal_error1_rank1_gb[i] = Fnorm(theta1_gas_rank1 - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_rank1_gb[i] = Fnorm(theta2_gas_rank1 - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  theta1_gas_center = tcrossprod(U0, V1) + tcrossprod(U1, A1)
  theta2_gas_center = trials * (tcrossprod(U0, V2) + tcrossprod(U2, A2))
  
  theta1_gas_center_rank1 = tcrossprod(U0_rank1, V1_rank1) + tcrossprod(U1_rank1, A1_rank1)
  theta2_gas_center_rank1 = trials * (tcrossprod(U0_rank1, V2_rank1) + tcrossprod(U2_rank1, A2_rank1))
  
  col1 = svd(theta1_gas_center)$u[,1:(r0+r1)]
  col2 = svd(theta2_gas_center)$u[,1:(r0+r2)]
  angle_result = angle_cal(col1, col2)
  
  col1_rank1 = svd(theta1_gas_center_rank1)$u[,1:(r0+r1)]
  col2_rank1 = svd(theta2_gas_center_rank1)$u[,1:(r0+r2)]
  angle_result_rank1 = angle_cal(col1_rank1, col2_rank1)
  
  # Get GAS canonical variables
  joint1_GAS = angle_result$principal_vector1[,1:r0]
  joint2_GAS = angle_result$principal_vector2[,1:r0]
  
  joint1_GAS_rank1 = angle_result_rank1$principal_vector1[,1:r0]
  joint2_GAS_rank1 = angle_result_rank1$principal_vector2[,1:r0]
  
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  GAS_joint_error1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS) - projection(joint1_true))
  GAS_joint_error2_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS) - projection(joint2_true))
  
  GAS_joint_error1_rank1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS_rank1) - projection(joint1_true))
  GAS_joint_error2_rank1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS_rank1) - projection(joint2_true))
}

###################### DCCA ##########################
for (i in 1:nrep){
  DCCA_signal_error1_gb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_gb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_gb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_gb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  DCCA_signal_error1_rank1_gb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_rank1_gb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_rank1_gb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_rank1_gb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
}

for (i in 1:nrep){
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  theta1_dcca_center = t(res_DCCA_list[[i]][[1]])
  theta2_dcca_center = t(res_DCCA_list[[i]][[2]])
  
  theta1_epca_dcca_center = t(res_EPCA_DCCA_list[[i]][[1]])
  theta2_epca_dcca_center = t(res_EPCA_DCCA_list[[i]][[2]])
  
  theta1_dcca_center_rank1 = t(res_DCCA_list_rank1[[i]][[1]])
  theta2_dcca_center_rank1 = t(res_DCCA_list_rank1[[i]][[2]])
  
  theta1_epca_dcca_center_rank1 = t(res_EPCA_DCCA_list_rank1[[i]][[1]])
  theta2_epca_dcca_center_rank1 = t(res_EPCA_DCCA_list_rank1[[i]][[2]])
  
  # DCCA joint
  col1_dcca = svd(theta1_dcca_center)$u[,1:(r0+r1)]
  col2_dcca = svd(theta2_dcca_center)$u[,1:(r0+r2)]
  angle_result_dcca = angle_cal(col1_dcca, col2_dcca)
  
  col1_dcca_rank1 = svd(theta1_dcca_center_rank1)$u[,1:(r0+r1)]
  col2_dcca_rank1 = svd(theta2_dcca_center_rank1)$u[,1:(r0+r2)]
  angle_result_dcca_rank1 = angle_cal(col1_dcca_rank1, col2_dcca_rank1)
  
  # EPCA-DCCA joint 
  col1_epca_dcca = svd(theta1_epca_dcca_center)$u[,1:(r0+r1)]
  col2_epca_dcca = svd(theta2_epca_dcca_center)$u[,1:(r0+r2)]
  angle_result_epca_dcca = angle_cal(col1_epca_dcca, col2_epca_dcca)
  
  col1_epca_dcca_rank1 = svd(theta1_epca_dcca_center_rank1)$u[,1:(r0+r1)]
  col2_epca_dcca_rank1 = svd(theta2_epca_dcca_center_rank1)$u[,1:(r0+r2)]
  angle_result_epca_dcca_rank1 = angle_cal(col1_epca_dcca_rank1, col2_epca_dcca_rank1)
  
  # Get DCCA canonical variables
  joint1_dcca = angle_result_dcca$principal_vector1[,1:r0]
  joint2_dcca = angle_result_dcca$principal_vector2[,1:r0]
  
  joint1_dcca_rank1 = angle_result_dcca_rank1$principal_vector1[,1:r0]
  joint2_dcca_rank1 = angle_result_dcca_rank1$principal_vector2[,1:r0]
  
  # Get EPCA-DCCA canonical variables
  joint1_epca_dcca = angle_result_epca_dcca$principal_vector1[,1:r0]
  joint2_epca_dcca = angle_result_epca_dcca$principal_vector2[,1:r0]
  
  joint1_epca_dcca_rank1 = angle_result_epca_dcca_rank1$principal_vector1[,1:r0]
  joint2_epca_dcca_rank1 = angle_result_epca_dcca_rank1$principal_vector2[,1:r0]
  
  #### Calculate the chordal distance between subspaces.
  DCCA_joint_error1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca) - projection(joint1_true))
  DCCA_joint_error2_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca) - projection(joint2_true))
  
  DCCA_joint_error1_rank1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca_rank1) - projection(joint1_true))
  DCCA_joint_error2_rank1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca_rank1) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca) - projection(joint1_true))
  EPCA_DCCA_joint_error2_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_rank1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca_rank1) - projection(joint1_true))
  EPCA_DCCA_joint_error2_rank1_gb[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca_rank1) - projection(joint2_true))
}

############################## Binomial - Binomial ##############################
source('Simulation/BB_Simulation_Trials100/GAS_BB_results.RData')
load('Simulation/BB_Simulation_Trials100/DCCA_result.RData')
load('Simulation/BB_Simulation_Trials100/ECCA_result.RData')
load('Simulation/BB_Simulation_Trials100/Data.RData')

# Initialize the signal error vectors.
ECCA_signal_error1_bb = rep(NA, nrep)
ECCA_signal_error2_bb = rep(NA, nrep)

DCCA_signal_error1_bb = rep(NA, nrep)
DCCA_signal_error2_bb = rep(NA, nrep)

EPCA_DCCA_signal_error1_bb = rep(NA, nrep)
EPCA_DCCA_signal_error2_bb = rep(NA, nrep)

GAS_signal_error1_bb = rep(NA, nrep)
GAS_signal_error2_bb = rep(NA, nrep)

DCCA_signal_error1_rank1_bb = rep(NA, nrep)
DCCA_signal_error2_rank1_bb = rep(NA, nrep)

EPCA_DCCA_signal_error1_rank1_bb = rep(NA, nrep)
EPCA_DCCA_signal_error2_rank1_bb = rep(NA, nrep)

GAS_signal_error1_rank1_bb = rep(NA, nrep)
GAS_signal_error2_rank1_bb = rep(NA, nrep)

###################### Compare the correlated signals ######################
# The idea is to get the first three canonical variables from estimated centered natural parameter matrix from GAS, DCCA, EPCA-DCCA 
# and compare it with the true 'correlated' signals.

# Initialize joint error vector.
ECCA_joint_error1_bb = rep(NA, nrep)
ECCA_joint_error2_bb = rep(NA, nrep)

DCCA_joint_error1_bb = rep(NA, nrep)
DCCA_joint_error2_bb = rep(NA, nrep)

EPCA_DCCA_joint_error1_bb = rep(NA, nrep)
EPCA_DCCA_joint_error2_bb = rep(NA, nrep)

GAS_joint_error1_bb = rep(NA, nrep)
GAS_joint_error2_bb = rep(NA, nrep)

# For joint rank 1 result:
DCCA_joint_error1_rank1_bb = rep(NA, nrep)
DCCA_joint_error2_rank1_bb = rep(NA, nrep)

EPCA_DCCA_joint_error1_rank1_bb = rep(NA, nrep)
EPCA_DCCA_joint_error2_rank1_bb = rep(NA, nrep)

GAS_joint_error1_rank1_bb = rep(NA, nrep)
GAS_joint_error2_rank1_bb = rep(NA, nrep)

###################### ECCA ##########################
for (i in 1:nrep){
  ECCA_signal_error1_bb[i] = Fnorm(output[[i]]$theta1_est - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  ECCA_signal_error2_bb[i] = Fnorm(output[[i]]$theta2_est - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  joint1_ecca = output[[i]]$U1_est
  joint2_ecca = output[[i]]$U2_est
  
  ECCA_joint_error1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_ecca) - projection(joint1_true))
  ECCA_joint_error2_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_ecca) - projection(joint2_true))
}

###################### GAS ##########################
for (i in 1:nrep){
  U0 = U0_GAS[,((i-1)*r0 + 1):(i*r0)]
  U1 = U1_GAS[,((i-1)*r1 + 1):(i*r1)]
  U2 = U2_GAS[,((i-1)*r2 + 1):(i*r2)]
  
  V1 = V0_GAS[1:p1,((i-1)*r0 + 1):(i*r0)]
  V2 = V0_GAS[(p1+1):(p1+p2),((i-1)*r0 + 1):(i*r0)]
  
  A1 = A1_GAS[,((i-1)*r1 + 1):(i*r1)]
  A2 = A2_GAS[,((i-1)*r2 + 1):(i*r2)]
  
  mu1 = Mu0_GAS[1:p1,i]
  mu2 = Mu0_GAS[(p1+1):(p1+p2),i]
  
  theta1_gas = trials * (tcrossprod(ones, mu1) + tcrossprod(U0, V1) + tcrossprod(U1, A1))
  theta2_gas = trials * (tcrossprod(ones, mu2) + tcrossprod(U0, V2) + tcrossprod(U2, A2))
  
  GAS_signal_error1_bb[i] = Fnorm(theta1_gas - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_bb[i] = Fnorm(theta2_gas - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  U0_rank1 = U0_GAS_rank1[,((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  U1_rank1 = U1_GAS_rank1[,((i-1)*r1_rank1 + 1):(i*r1_rank1)]
  U2_rank1 = U2_GAS_rank1[,((i-1)*r2_rank1 + 1):(i*r2_rank1)]
  
  V1_rank1 = V0_GAS_rank1[1:p1,((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  V2_rank1 = V0_GAS_rank1[(p1+1):(p1+p2),((i-1)*r0_rank1 + 1):(i*r0_rank1)]
  
  A1_rank1 = A1_GAS_rank1[,((i-1)*r1_rank1 + 1):(i*r1_rank1)]
  A2_rank1 = A2_GAS_rank1[,((i-1)*r2_rank1 + 1):(i*r2_rank1)]
  
  mu1_rank1 = Mu0_GAS_rank1[1:p1,i]
  mu2_rank1 = Mu0_GAS_rank1[(p1+1):(p1+p2),i]
  
  theta1_gas_rank1 = trials * (tcrossprod(ones, mu1_rank1) + tcrossprod(U0_rank1, V1_rank1) + tcrossprod(U1_rank1, A1_rank1))
  theta2_gas_rank1 = trials * (tcrossprod(ones, mu2_rank1) + tcrossprod(U0_rank1, V2_rank1) + tcrossprod(U2_rank1, A2_rank1))
  
  GAS_signal_error1_rank1_bb[i] = Fnorm(theta1_gas_rank1 - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_rank1_bb[i] = Fnorm(theta2_gas_rank1 - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  theta1_gas_center = trials * (tcrossprod(U0, V1) + tcrossprod(U1, A1))
  theta2_gas_center = trials * (tcrossprod(U0, V2) + tcrossprod(U2, A2))
  
  theta1_gas_center_rank1 = trials * (tcrossprod(U0_rank1, V1_rank1) + tcrossprod(U1_rank1, A1_rank1))
  theta2_gas_center_rank1 = trials * (tcrossprod(U0_rank1, V2_rank1) + tcrossprod(U2_rank1, A2_rank1))
  
  col1 = svd(theta1_gas_center)$u[,1:(r0+r1)]
  col2 = svd(theta2_gas_center)$u[,1:(r0+r2)]
  angle_result = angle_cal(col1, col2)
  
  col1_rank1 = svd(theta1_gas_center_rank1)$u[,1:(r0+r1)]
  col2_rank1 = svd(theta2_gas_center_rank1)$u[,1:(r0+r2)]
  angle_result_rank1 = angle_cal(col1_rank1, col2_rank1)
  
  # Get GAS canonical variables
  joint1_GAS = angle_result$principal_vector1[,1:r0]
  joint2_GAS = angle_result$principal_vector2[,1:r0]
  
  joint1_GAS_rank1 = angle_result_rank1$principal_vector1[,1:r0]
  joint2_GAS_rank1 = angle_result_rank1$principal_vector2[,1:r0]
  
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  GAS_joint_error1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS) - projection(joint1_true))
  GAS_joint_error2_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS) - projection(joint2_true))
  
  GAS_joint_error1_rank1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS_rank1) - projection(joint1_true))
  GAS_joint_error2_rank1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS_rank1) - projection(joint2_true))
}

###################### DCCA ##########################
for (i in 1:nrep){
  DCCA_signal_error1_bb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_bb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_bb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_bb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  DCCA_signal_error1_rank1_bb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_rank1_bb[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_rank1_bb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_rank1_bb[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
}

for (i in 1:nrep){
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  theta1_dcca_center = t(res_DCCA_list[[i]][[1]])
  theta2_dcca_center = t(res_DCCA_list[[i]][[2]])
  
  theta1_epca_dcca_center = t(res_EPCA_DCCA_list[[i]][[1]])
  theta2_epca_dcca_center = t(res_EPCA_DCCA_list[[i]][[2]])
  
  theta1_dcca_center_rank1 = t(res_DCCA_list_rank1[[i]][[1]])
  theta2_dcca_center_rank1 = t(res_DCCA_list_rank1[[i]][[2]])
  
  theta1_epca_dcca_center_rank1 = t(res_EPCA_DCCA_list_rank1[[i]][[1]])
  theta2_epca_dcca_center_rank1 = t(res_EPCA_DCCA_list_rank1[[i]][[2]])
  
  # DCCA joint
  col1_dcca = svd(theta1_dcca_center)$u[,1:(r0+r1)]
  col2_dcca = svd(theta2_dcca_center)$u[,1:(r0+r2)]
  angle_result_dcca = angle_cal(col1_dcca, col2_dcca)
  
  col1_dcca_rank1 = svd(theta1_dcca_center_rank1)$u[,1:(r0+r1)]
  col2_dcca_rank1 = svd(theta2_dcca_center_rank1)$u[,1:(r0+r2)]
  angle_result_dcca_rank1 = angle_cal(col1_dcca_rank1, col2_dcca_rank1)
  
  # EPCA-DCCA joint 
  col1_epca_dcca = svd(theta1_epca_dcca_center)$u[,1:(r0+r1)]
  col2_epca_dcca = svd(theta2_epca_dcca_center)$u[,1:(r0+r2)]
  angle_result_epca_dcca = angle_cal(col1_epca_dcca, col2_epca_dcca)
  
  col1_epca_dcca_rank1 = svd(theta1_epca_dcca_center_rank1)$u[,1:(r0+r1)]
  col2_epca_dcca_rank1 = svd(theta2_epca_dcca_center_rank1)$u[,1:(r0+r2)]
  angle_result_epca_dcca_rank1 = angle_cal(col1_epca_dcca_rank1, col2_epca_dcca_rank1)
  
  # Get DCCA canonical variables
  joint1_dcca = angle_result_dcca$principal_vector1[,1:r0]
  joint2_dcca = angle_result_dcca$principal_vector2[,1:r0]
  
  joint1_dcca_rank1 = angle_result_dcca_rank1$principal_vector1[,1:r0]
  joint2_dcca_rank1 = angle_result_dcca_rank1$principal_vector2[,1:r0]
  
  # Get EPCA-DCCA canonical variables
  joint1_epca_dcca = angle_result_epca_dcca$principal_vector1[,1:r0]
  joint2_epca_dcca = angle_result_epca_dcca$principal_vector2[,1:r0]
  
  joint1_epca_dcca_rank1 = angle_result_epca_dcca_rank1$principal_vector1[,1:r0]
  joint2_epca_dcca_rank1 = angle_result_epca_dcca_rank1$principal_vector2[,1:r0]
  
  #### Calculate the chordal distance between subspaces.
  DCCA_joint_error1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca) - projection(joint1_true))
  DCCA_joint_error2_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca) - projection(joint2_true))
  
  DCCA_joint_error1_rank1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca_rank1) - projection(joint1_true))
  DCCA_joint_error2_rank1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca_rank1) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca) - projection(joint1_true))
  EPCA_DCCA_joint_error2_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_rank1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca_rank1) - projection(joint1_true))
  EPCA_DCCA_joint_error2_rank1_bb[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca_rank1) - projection(joint2_true))
}

# ggplot for the signal error:
library(ggplot2)

method = c(rep('ECCA', 6*nrep), rep('DCCA', 6*nrep), rep('EPCA-DCCA', 6*nrep), rep('GAS-rank3', 6*nrep), rep('GAS-rank1', 6*nrep))
setting = rep(c(rep('G-G 1st signal', nrep), rep('G-G 2nd signal', nrep),
                rep('G-B 1st signal', nrep), rep('G-B 2nd signal', nrep),
                rep('B-B 1st signal', nrep), rep('B-B 2nd signal', nrep)),5)
y_error = c(ECCA_signal_error1_gg, ECCA_signal_error2_gg,
            ECCA_signal_error1_gb, ECCA_signal_error2_gb,
            ECCA_signal_error1_bb, ECCA_signal_error2_bb,
            
            DCCA_signal_error1_gg, DCCA_signal_error2_gg, 
            DCCA_signal_error1_gb, DCCA_signal_error2_gb,
            DCCA_signal_error1_bb, DCCA_signal_error2_bb,
            
            EPCA_DCCA_signal_error1_gg, EPCA_DCCA_signal_error2_gg,
            EPCA_DCCA_signal_error1_gb, EPCA_DCCA_signal_error2_gb,
            EPCA_DCCA_signal_error1_bb, EPCA_DCCA_signal_error2_bb,
            
            GAS_signal_error1_gg, GAS_signal_error2_gg, 
            GAS_signal_error1_gb, GAS_signal_error2_gb,
            GAS_signal_error1_bb, GAS_signal_error2_bb,
            
            GAS_signal_error1_rank1_gg, GAS_signal_error2_rank1_gg,
            GAS_signal_error1_rank1_gb, GAS_signal_error2_rank1_gb,
            GAS_signal_error1_rank1_bb, GAS_signal_error2_rank1_bb)

gg_mat = data.frame(Error = y_error, Method = method, setting = setting)

gg_mat$setting <- factor(gg_mat$setting,
                         c('G-G 1st signal', 'G-B 1st signal', 'B-B 1st signal', 'G-G 2nd signal', 'G-B 2nd signal', 'B-B 2nd signal'))

gg <- ggplot(gg_mat, aes(x=Method, y=Error))
gg <- gg + geom_boxplot(aes(color=Method))
gg <- gg + facet_wrap(~setting)
gg <- gg + theme_bw()
gg <- gg + ylab("Relative Error")
gg <- gg + theme(strip.background=element_rect(fill="black"))
gg <- gg + theme(strip.text=element_text(color="white", face="bold",size = 25))
gg <- gg + theme(text=element_text(size = 25))
gg <- gg + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
gg <- gg + scale_colour_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))
gg

# ggplot for the joint error:
y_joint = c(ECCA_joint_error1_gg, ECCA_joint_error2_gg,
            ECCA_joint_error1_gb, ECCA_joint_error2_gb,
            ECCA_joint_error1_bb, ECCA_joint_error2_bb,
            
            DCCA_joint_error1_gg, DCCA_joint_error2_gg, 
            DCCA_joint_error1_gb, DCCA_joint_error2_gb,
            DCCA_joint_error1_bb, DCCA_joint_error2_bb,
            
            EPCA_DCCA_joint_error1_gg, EPCA_DCCA_joint_error2_gg,
            EPCA_DCCA_joint_error1_gb, EPCA_DCCA_joint_error2_gb,
            EPCA_DCCA_joint_error1_bb, EPCA_DCCA_joint_error2_bb,
            
            GAS_joint_error1_gg, GAS_joint_error2_gg, 
            GAS_joint_error1_gb, GAS_joint_error2_gb,
            GAS_joint_error1_bb, GAS_joint_error2_bb,
            
            GAS_joint_error1_rank1_gg, GAS_joint_error2_rank1_gg,
            GAS_joint_error1_rank1_gb, GAS_joint_error2_rank1_gb,
            GAS_joint_error1_rank1_bb, GAS_joint_error2_rank1_bb)

setting_joint = rep(c(rep('G-G 1st joint', nrep), rep('G-G 2nd joint', nrep),
                      rep('G-B 1st joint', nrep), rep('G-B 2nd joint', nrep),
                      rep('B-B 1st joint', nrep), rep('B-B 2nd joint', nrep)),5)

gg_joint_mat = data.frame(Error = y_joint, Method = method, setting_joint = setting_joint)

gg_joint_mat$setting_joint <- factor(gg_joint_mat$setting_joint,
                         c('G-G 1st joint', 'G-B 1st joint', 'B-B 1st joint', 'G-G 2nd joint', 'G-B 2nd joint', 'B-B 2nd joint'))

gg_joint <- ggplot(gg_joint_mat, aes(x=Method, y=Error))
gg_joint <- gg_joint + geom_boxplot(aes(color=Method))
gg_joint <- gg_joint + facet_wrap(~setting_joint)
gg_joint <- gg_joint + theme_bw()
gg_joint <- gg_joint + ylab("Chordal Distance")
gg_joint <- gg_joint + theme(strip.background=element_rect(fill="black"))
gg_joint <- gg_joint + theme(strip.text=element_text(color="white", face="bold",size = 25))
gg_joint <- gg_joint + theme(text=element_text(size = 25))
gg_joint <- gg_joint + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
gg_joint <- gg_joint + scale_colour_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))
gg_joint

fig.path = "Simulation/Figures/"
pdf(file = paste(fig.path,"Whole_Estimation_Error.pdf",sep=""), width = 16, height = 9)
print(gg)
dev.off()

fig.path = "Simulation/Figures/"
pdf(file = paste(fig.path,"Whole_Joint_Error.pdf",sep=""), width = 16, height = 9)
print(gg_joint)
dev.off()
