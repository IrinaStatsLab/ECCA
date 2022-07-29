rm(list = ls())

source('Function/exp_family.R')
source('Simulation/GG_Simulation_SNR5/GAS_GG_results.RData')
load('Simulation/GG_Simulation_SNR5/DCCA_result.RData')
load('Simulation/GG_Simulation_SNR5/ECCA_result.RData')
load('Simulation/GG_Simulation_SNR5/Data.RData')

# Set parameters
n = 50
p1 = 30
p2 = 20
p = c(p1, p2)
family1 = 'gaussian'
family2 = 'gaussian'
lambdas = c(1, 0.9, 0.7) # covariance between joint structures

# Set ranks (joint and individuals)
r0 = 3
r1 = 4
r2 = 3

r0_rank1 = 1
r1_rank1 = 6
r2_rank1 = 5

nrep = 100
set.seed(37)

# Initialize the signal error vectors.
ECCA_signal_error1 = rep(NA, nrep)
ECCA_signal_error2 = rep(NA, nrep)

DCCA_signal_error1 = rep(NA, nrep)
DCCA_signal_error2 = rep(NA, nrep)

EPCA_DCCA_signal_error1 = rep(NA, nrep)
EPCA_DCCA_signal_error2 = rep(NA, nrep)

GAS_signal_error1 = rep(NA, nrep)
GAS_signal_error2 = rep(NA, nrep)

DCCA_signal_error1_rank1 = rep(NA, nrep)
DCCA_signal_error2_rank1 = rep(NA, nrep)

EPCA_DCCA_signal_error1_rank1 = rep(NA, nrep)
EPCA_DCCA_signal_error2_rank1 = rep(NA, nrep)

GAS_signal_error1_rank1 = rep(NA, nrep)
GAS_signal_error2_rank1 = rep(NA, nrep)

###################### Compare the correlated signals ######################
# The idea is to get the first three canonical variables from estimated centered natural parameter matrix from GAS, DCCA, EPCA-DCCA 
# and compare it with the true 'correlated' signals.

# Initialize joint error vector.
ECCA_joint_error1 = rep(NA, nrep)
ECCA_joint_error2 = rep(NA, nrep)

DCCA_joint_error1 = rep(NA, nrep)
DCCA_joint_error2 = rep(NA, nrep)

EPCA_DCCA_joint_error1 = rep(NA, nrep)
EPCA_DCCA_joint_error2 = rep(NA, nrep)

GAS_joint_error1 = rep(NA, nrep)
GAS_joint_error2 = rep(NA, nrep)

# For joint rank 1 result:
DCCA_joint_error1_rank1 = rep(NA, nrep)
DCCA_joint_error2_rank1 = rep(NA, nrep)

EPCA_DCCA_joint_error1_rank1 = rep(NA, nrep)
EPCA_DCCA_joint_error2_rank1 = rep(NA, nrep)

GAS_joint_error1_rank1 = rep(NA, nrep)
GAS_joint_error2_rank1 = rep(NA, nrep)

# Investigate the estimated first three cc:
ECCA_cc_mat = matrix(rep(0, 3*nrep), nrow = nrep)
DCCA_cc_mat = matrix(rep(0, 3*nrep), nrow = nrep)
EPCA_DCCA_cc_mat = matrix(rep(0, 3*nrep), nrow = nrep)
GAS_cc_mat = matrix(rep(0, 3*nrep), nrow = nrep)
GAS_cc_mat_rank1 = matrix(rep(0, 3*nrep), nrow = nrep)

ones = rep(1, n)
###################### ECCA ##########################
for (i in 1:nrep){
  ECCA_signal_error1[i] = Fnorm(theta1_est[[i]] - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  ECCA_signal_error2[i] = Fnorm(theta2_est[[i]] - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  ######### True signals #########
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  joint1_ecca = U1_est[[i]]
  joint2_ecca = U2_est[[i]]
  
  ECCA_joint_error1[i] = 1/sqrt(2)*Fnorm(projection(joint1_ecca) - projection(joint1_true))
  ECCA_joint_error2[i] = 1/sqrt(2)*Fnorm(projection(joint2_ecca) - projection(joint2_true))
  
  temp_angle = angle_cal(joint1_ecca, joint2_ecca)$angle
  ECCA_cc_mat[i,] = cos(temp_angle)
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
  
  GAS_signal_error1[i] = Fnorm(theta1_gas - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2[i] = Fnorm(theta2_gas - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
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
  
  GAS_signal_error1_rank1[i] = Fnorm(theta1_gas_rank1 - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  GAS_signal_error2_rank1[i] = Fnorm(theta2_gas_rank1 - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
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
  # print(cos(angle_result$angle))
  # print(cos(angle_result_rank1$angle))
  
  # Get GAS canonical variables
  joint1_GAS = angle_result$principal_vector1[,1:r0]
  joint2_GAS = angle_result$principal_vector2[,1:r0]
  
  joint1_GAS_rank1 = angle_result_rank1$principal_vector1[,1:r0]
  joint2_GAS_rank1 = angle_result_rank1$principal_vector2[,1:r0]
  
  joint1_true = U1_true[[i]]
  joint2_true = U2_true[[i]]
  
  GAS_joint_error1[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS) - projection(joint1_true))
  GAS_joint_error2[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS) - projection(joint2_true))
  
  GAS_joint_error1_rank1[i] = 1/sqrt(2)*Fnorm(projection(joint1_GAS_rank1) - projection(joint1_true))
  GAS_joint_error2_rank1[i] = 1/sqrt(2)*Fnorm(projection(joint2_GAS_rank1) - projection(joint2_true))
  
  temp_angle = angle_result$angle
  GAS_cc_mat[i,] = cos(temp_angle)[1:3]
  
  temp_angle_rank1 = angle_result_rank1$angle
  GAS_cc_mat_rank1[i,] = cos(temp_angle_rank1)[1:3]
}

###################### DCCA ##########################
for (i in 1:nrep){
  DCCA_signal_error1[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  DCCA_signal_error1_rank1[i] = Fnorm(tcrossprod(ones, res_DCCA_mu1_list[[i]]) + t(res_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  DCCA_signal_error2_rank1[i] = Fnorm(tcrossprod(ones, res_DCCA_mu2_list[[i]]) + t(res_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
  
  EPCA_DCCA_signal_error1_rank1[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu1_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[1]]) - theta1_true[[i]])^2/Fnorm(theta1_true[[i]])^2
  EPCA_DCCA_signal_error2_rank1[i] = Fnorm(tcrossprod(ones, res_EPCA_DCCA_mu2_list[[i]]) + t(res_EPCA_DCCA_list_rank1[[i]][[2]]) - theta2_true[[i]])^2/Fnorm(theta2_true[[i]])^2
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
  DCCA_joint_error1[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca) - projection(joint1_true))
  DCCA_joint_error2[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca) - projection(joint2_true))
  
  DCCA_joint_error1_rank1[i] = 1/sqrt(2)*Fnorm(projection(joint1_dcca_rank1) - projection(joint1_true))
  DCCA_joint_error2_rank1[i] = 1/sqrt(2)*Fnorm(projection(joint2_dcca_rank1) - projection(joint2_true))

  EPCA_DCCA_joint_error1[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca) - projection(joint1_true))
  EPCA_DCCA_joint_error2[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca) - projection(joint2_true))
  
  EPCA_DCCA_joint_error1_rank1[i] = 1/sqrt(2)*Fnorm(projection(joint1_epca_dcca_rank1) - projection(joint1_true))
  EPCA_DCCA_joint_error2_rank1[i] = 1/sqrt(2)*Fnorm(projection(joint2_epca_dcca_rank1) - projection(joint2_true))
  
  temp_angle_dcca = angle_result_dcca$angle
  DCCA_cc_mat[i,] = cos(temp_angle_dcca)[1:3]
  
  temp_angle_epca_dcca = angle_result_epca_dcca$angle
  EPCA_DCCA_cc_mat[i,] = cos(temp_angle_epca_dcca)[1:3]
}

boxplot(ECCA_signal_error1, DCCA_signal_error1, EPCA_DCCA_signal_error1, GAS_signal_error1, DCCA_signal_error1_rank1, EPCA_DCCA_signal_error1_rank1, GAS_signal_error1_rank1)
boxplot(ECCA_signal_error2, DCCA_signal_error2, EPCA_DCCA_signal_error2, GAS_signal_error2, DCCA_signal_error2_rank1, EPCA_DCCA_signal_error2_rank1, GAS_signal_error2_rank1)

boxplot(ECCA_joint_error1, DCCA_joint_error1, EPCA_DCCA_joint_error1, GAS_joint_error1, DCCA_joint_error1_rank1, EPCA_DCCA_joint_error1_rank1, GAS_joint_error1_rank1)
boxplot(ECCA_joint_error2, DCCA_joint_error2, EPCA_DCCA_joint_error2, GAS_joint_error2, DCCA_joint_error2_rank1, EPCA_DCCA_joint_error2_rank1, GAS_joint_error2_rank1)

# Let's check the average cc for each method:
colMeans(ECCA_cc_mat)
colMeans(DCCA_cc_mat)
colMeans(EPCA_DCCA_cc_mat)
colMeans(GAS_cc_mat)
colMeans(GAS_cc_mat_rank1)
boxplot(ECCA_cc_mat[,1])

# Get the standard deviation:
library(GMCM)
GMCM:::colSds(ECCA_cc_mat)
# sd(ECCA_cc_mat[,1]) # sanity check
GMCM:::colSds(DCCA_cc_mat)
GMCM:::colSds(EPCA_DCCA_cc_mat)
GMCM:::colSds(GAS_cc_mat)
GMCM:::colSds(GAS_cc_mat_rank1)

# ggplot for the signal error:
library(ggplot2)

method = c(rep('ECCA', 2*nrep), rep('DCCA', 2*nrep), rep('EPCA-DCCA', 2*nrep), rep('GAS-rank3', 2*nrep), rep('GAS-rank1', 2*nrep))
matrix = rep(c(rep('1st Signal', nrep), rep('2nd Signal', nrep)),5)
y_error = c(ECCA_signal_error1, ECCA_signal_error2,
            DCCA_signal_error1, DCCA_signal_error2, 
            EPCA_DCCA_signal_error1, EPCA_DCCA_signal_error2,
            GAS_signal_error1, GAS_signal_error2, 
            GAS_signal_error1_rank1, GAS_signal_error2_rank1)

gg_mat = data.frame(Error = y_error, Method = method, Section = matrix)

gg <- ggplot(gg_mat, aes(x=Method, y=Error))
gg <- gg + geom_boxplot(aes(color=Method))
gg <- gg + facet_wrap(~Section)
gg <- gg + theme_bw()
gg <- gg + ylab("Relative Error")
gg <- gg + theme(strip.background=element_rect(fill="black"))
gg <- gg + theme(strip.text=element_text(color="white", face="bold",size = 25))
gg <- gg + theme(text=element_text(size = 25))
gg <- gg + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
gg <- gg + scale_colour_manual(values = c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))
gg

# ggplot for the joint error:
matrix_joint = rep(c(rep('1st Joint', nrep), rep('2nd Joint', nrep)),5)
y_error_joint = c(ECCA_joint_error1, ECCA_joint_error2,
                  DCCA_joint_error1, DCCA_joint_error2, 
                  EPCA_DCCA_joint_error1, EPCA_DCCA_joint_error2, 
                  GAS_joint_error1, GAS_joint_error2, 
                  GAS_joint_error1_rank1, GAS_joint_error2_rank1)

gg_mat_joint = data.frame(Error = y_error_joint, Method = method, Section = matrix_joint)

gg_joint <- ggplot(gg_mat_joint, aes(x=Method, y=Error))
gg_joint <- gg_joint + geom_boxplot(aes(color=Method))
gg_joint <- gg_joint + facet_wrap(~Section)
gg_joint <- gg_joint + theme_bw()
gg_joint <- gg_joint + ylab("Chordal Distance")
gg_joint <- gg_joint + theme(strip.background=element_rect(fill="black"))
gg_joint <- gg_joint + theme(strip.text=element_text(color="white", face="bold",size = 25))
gg_joint <- gg_joint + theme(text=element_text(size = 25))
gg_joint <- gg_joint + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
gg_joint <- gg_joint + scale_colour_manual(values =  c('#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00'))
gg_joint

fig.path = "Simulation/GG_Simulation_SNR5/Figures/"
pdf(file = paste(fig.path,"Whole_Estimation_Error_GG_SNR5_Draft.pdf",sep=""), width = 16, height = 9)
print(gg)
dev.off()

fig.path = "Simulation/GG_Simulation_SNR5/Figures/"
pdf(file = paste(fig.path,"Whole_Joint_Error_GG_SNR5_subspace_Draft.pdf",sep=""), width = 16, height = 9)
print(gg_joint)
dev.off()