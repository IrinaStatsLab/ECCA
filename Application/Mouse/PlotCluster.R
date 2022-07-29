# This file is used to generate the graphs shown in the manuscript

rm(list = ls())
load("Application/Mouse/ECCA_result.RData")
library(mixOmics)
data(nutrimouse)
library(ggplot2)
library(tidyverse)

# total ranks r1 = 3 and r2 = 4 
# joint rank r0 = 2

out100 = svd(eccaMat100$Z2 %*% t(eccaMat100$A2))$u[, 1:2]
out100_joint = svd(cbind(eccaMat100$U1, eccaMat100$U2))$u[,1:2]

newData100 = data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype, indiv1 = out100[, 1], indiv2 = out100[, 2])

gg100 = newData100 %>%
  ggplot(aes(x = indiv1, y = indiv2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(gg100)

fig.path = "Application/Mouse/Figures/"
pdf(file = paste(fig.path,"ECCA_Ind_Separation_Mouse_not_Scaling_trials100.pdf",sep=""), width = 10, height = 8)
print(gg100)
dev.off()

newData100_joint = data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype, joint1 = out100_joint[, 1], joint2 = out100_joint[, 2])

gg100_joint = newData100_joint %>%
  ggplot(aes(x = joint1, y = joint2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(gg100_joint)

pdf(file = paste(fig.path,"ECCA_Joint_Separation_Mouse_not_Scaling_trials100.pdf",sep=""), width = 10, height = 8)
print(gg100_joint)
dev.off()

#################### GAS ##########################
source("Application/Mouse/GAS_result.RData")
# Source the angle_cal function
source('Function/exp_family.R')
# We report the result of setup 1, which is the same rank setup as ECCA here
U0 = U0_set1
U1 = U1_set1
U2 = U2_set1

r0 = 2
r1 = 1
r2 = 2
p1 = dim(nutrimouse$gene)[2] # 120
p2 = dim(nutrimouse$lipid)[2] # 21
trials = 100
V1 = V0_set1[1:p1,1:r0]
V2 = V0_set1[(p1+1):(p1+p2),1:r0]

A1 = A1_set1
A2 = A2_set1

# Calculate the angles between the column spaces of individual score matrices:
# individual ranks: r1 = 1, r2 = 2, so the result will only contain one angle
angle_result = angle_cal(U1, U2)

angle_result$angle # 0.91 radian, around 52 degrees

#### Plot the similar graph as ECCA ####
# Individual
out_GAS = svd(U2 %*% t(A2))$u[, 1:2]
# Joint
out_jointGAS = svd(U0)$u[,1:2]

newData_GAS = data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype, indiv1 = out_GAS[, 1], indiv2 = out_GAS[, 2])

newData_jointGAS = data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype, joint1 = out_jointGAS[, 1], joint2 = out_jointGAS[, 2])

gg_GAS = newData_GAS %>%
  ggplot(aes(x = indiv1, y = indiv2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(gg_GAS)

gg_jointGAS = newData_jointGAS %>%
  ggplot(aes(x = joint1, y = joint2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(gg_jointGAS)

fig.path = "Application/Mouse/Figures/"
pdf(file = paste(fig.path,"GAS_Ind_Separation_Mouse.pdf",sep = ""), width = 10, height = 8)
print(gg_GAS)
dev.off()

fig.path = "Application/Mouse/Figures/"
pdf(file = paste(fig.path,"GAS_Joint_Separation_Mouse.pdf",sep = ""), width = 10, height = 8)
print(gg_jointGAS)
dev.off()