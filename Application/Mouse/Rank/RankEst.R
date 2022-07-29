# This file is used to determine the ranks. 
# Result: r0 = 2 (joint rank), r1 = 1 (individual rank 1), r2 = 2 (individual rank 2) based on BCV and principal angles.

rm(list = ls())
# https://mixomicsteam.github.io/Bookdown/ for the tutorial for mixOmics
library(mixOmics)
?nutrimouse
data(nutrimouse)
require(pracma)
library('R.matlab')
library('matlabr')

x1 = as.matrix(nutrimouse$gene)
# x1 = scale(x1)
x2 = as.matrix(nutrimouse$lipid/100)
# We need to use some tricks to make sure our algorithm works on this data, since there are zeros.
# Wu use the trick mentioned in Dr. Longnecker's book:
# if x = 0 then set it to be 3/8/(100+3/4)
for (i in 1:dim(x2)[1]){
  for (j in 1:dim(x2)[2]){
    if (x2[i,j] == 0){
      x2[i,j] = 3/8/(100+3/4)
    }
  }
}
# family1 = 'gaussian'
# family2 = 'binomial'

# Save as the mat file
writeMat('Application/Mouse/Rank/Mouse.mat', x1 = x1, x2 = x2)
# Use the rank estimation procedure of GAS to estimate the ranks of Its and Timer.

# Takes 10-20 minutes to run.
run_matlab_script('Application/Mouse/Rank/RankEst.m')
source('Application/Mouse/Rank/MouseRank.RData')
r1 # 3
r2 # 7 might be too large. Let's investigate.
plot(apply(temp_allCVscore1,2,median)[1:10], xlab = 'Rank', ylab = 'Median CV score for x1')
plot(apply(temp_allCVscore2,2,median)[1:10], xlab = 'Rank', ylab = 'Median CV score for x2') # r_2 >= 4.
plot(apply(temp_allCVscore2,2,mean)[1:10], xlab = 'Rank', ylab = 'Mean CV score for x2')
plot(apply(temp_allCVscore2,2,mean)[4:10], xlab = 'Rank', ylab = 'Median CV score for x2')
# r2 ranges from 4 to 7. It seems that rank 4 claim is close to the error of rank 7 claim.
# In summary, we use total rank 1 to be 3, and total rank 2 to be 4.

# Use exponential PCA to determine joint rank:
library(generalizedPCA)
source('Function/exp_family.R')
# Total ranks are 3 and 4 respectively.
# Get the low-rank natural parameter from exponential PCA:
para1 = generalizedPCA(x1, k = 3, family = 'gaussian')
para2 = generalizedPCA(x2, k = 4, family = 'binomial')

ones = rep(1, dim(x1)[1])
theta1_epca = (x1 - tcrossprod(ones, para1$mu)) %*% tcrossprod(para1$U)
theta2_epca = (x2 - tcrossprod(ones, para2$mu)) %*% tcrossprod(para2$U)

theta1_epca_c = as.matrix(svd(theta1_epca)$u[,1:3])
theta2_epca_c = as.matrix(svd(theta2_epca)$u[,1:4])

# Calculate the principal angles.
angle_cal(theta1_epca_c, theta2_epca_c)$angle * 180/pi
# 74 degree is close to 90. Set joint rank as 2.