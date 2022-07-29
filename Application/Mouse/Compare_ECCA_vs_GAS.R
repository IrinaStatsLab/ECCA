# Compare results of GAS with results of ECCA
# SWISS score - way to evaluate cluster consistency
# used in JIVE Lock et al., 2011 paper, originally described in Cabanski et al., 2010

# Function to calculate SWISS score
# X - n times p matrix of measurements
# class - a priori cluster/subtype information for n samples, integer from 1 to C
###################################################################################
# On JIVE GBM data, the function results agree with Table 2 results in Lock et al., 2011 as long as unclassified samples (subtypes == 0) are removed
calculateSwiss <- function(X, class){
  # Calculate total sum of squares
  ################################
  # Subtract column means 
  Xs <- scale(X, scale = F)
  SST <- sum(Xs^2)
  # Calculate within-class sum of squares
  TWISS <- 0
  for (i in 1:max(class)){
    TWISS <- TWISS + sum((X[class == i, , drop = FALSE] - matrix(colMeans(X[class == i, , drop = FALSE]), sum(class == i), ncol(X), byrow = T))^2)
  }
  return(TWISS/SST)
}


# Load results
library(R.matlab)
GASresults = readMat('Application/Mouse/FitGAS/Mice_GAS_output.mat')
load("Application/Mouse/ECCA_result.RData")

# Joint ECCA versus Joint GAS
source("Function/exp_family.R")
# ECCA joint correlations 0.87 and 0.65 
round(crossprod(eccaMat100$U1, eccaMat100$U2), 2)  
# versus GAS being 1 and 1 (just U0.set1)

# ECCA individual correlations are 0 and 0
round(crossprod(eccaMat100$Z1, eccaMat100$Z2), 2)  

# GAS individual correlations are -0.16 and 0.59
round(cor(GASresults$U1.set1, GASresults$U2.set1), 2)

# Load data
library(mixOmics)
data(nutrimouse)

# Plots
library(tidyverse)

# Rotate ECCA individual within diet for identifiability
ECCAindiv_lipids = svd(eccaMat100$Z2 %*% t(eccaMat100$A2))$u[, 1:2]
# Combine two joint for ECCA
ECCAjoint = svd(cbind(eccaMat100$U1, eccaMat100$U2))$u[,1:2]
# Standardize GAS joint
GASjoint = svd(GASresults$U0.set1 %*% t(GASresults$V0.set1))$u[, 1:2]
# Rotate GAS individual as for ECCA
GASindiv_lipids = svd(GASresults$U2.set1 %*% t(GASresults$A2.set1))$u[, 1:2]

# Calculate SWISS scores on joint structures
# For genotype
calculateSwiss(GASjoint, as.numeric(nutrimouse$genotype)) #0.626
calculateSwiss(ECCAjoint, as.numeric(nutrimouse$genotype)) #0.571
# Looking at one at a time
calculateSwiss(GASjoint[,1, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.41
calculateSwiss(GASjoint[,2, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.84
calculateSwiss(ECCAjoint[,1, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.145
calculateSwiss(ECCAjoint[,2, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.998

# For diet
calculateSwiss(GASjoint, as.numeric(nutrimouse$diet)) #0.854
calculateSwiss(ECCAjoint, as.numeric(nutrimouse$diet)) #0.587
# Looking at one at a time
calculateSwiss(GASjoint[,1, drop = FALSE], as.numeric(nutrimouse$diet)) #0.85
calculateSwiss(GASjoint[,2, drop = FALSE], as.numeric(nutrimouse$diet)) #0.86
calculateSwiss(ECCAjoint[,1, drop = FALSE], as.numeric(nutrimouse$diet)) #0.94
calculateSwiss(ECCAjoint[,2, drop = FALSE], as.numeric(nutrimouse$diet)) #0.23

# Calculate SWISS scores on individual structures
# For genotype
calculateSwiss(GASindiv_lipids, as.numeric(nutrimouse$genotype)) #0.95
calculateSwiss(ECCAindiv_lipids, as.numeric(nutrimouse$genotype)) #0.98
# Looking at one at a time
calculateSwiss(GASindiv_lipids[,1, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.99
calculateSwiss(GASindiv_lipids[,2, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.91
calculateSwiss(ECCAindiv_lipids[,1, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.99
calculateSwiss(ECCAindiv_lipids[,2, drop = FALSE], as.numeric(nutrimouse$genotype)) #0.96

# For diet
calculateSwiss(GASindiv_lipids, as.numeric(nutrimouse$diet)) #0.152
calculateSwiss(ECCAindiv_lipids, as.numeric(nutrimouse$diet)) #0.154
# Looking at one at a time
calculateSwiss(GASindiv_lipids[,1, drop = FALSE], as.numeric(nutrimouse$diet)) #0.11
calculateSwiss(GASindiv_lipids[,2, drop = FALSE], as.numeric(nutrimouse$diet)) #0.19
calculateSwiss(ECCAindiv_lipids[,1, drop = FALSE], as.numeric(nutrimouse$diet)) #0.06
calculateSwiss(ECCAindiv_lipids[,2, drop = FALSE], as.numeric(nutrimouse$diet)) #0.25

# Combine into one data frame for plotting (see PlotCluster.R for saved plots)
plotData = data.frame(diet = nutrimouse$diet, genotype = nutrimouse$genotype,
                      GASjoint1 = GASjoint[, 1], GASjoint2 = GASjoint[, 2],
                      ECCAjoint1 = ECCAjoint[, 1], ECCAjoint2 = ECCAjoint[, 2],
                      ECCAindiv1 = ECCAindiv_lipids[, 1], ECCAindiv2 = ECCAindiv_lipids[, 2],
                      GASindiv1 = GASindiv_lipids[, 1], GASindiv2 = GASindiv_lipids[, 2])

p1 = plotData %>%
  ggplot(aes(x = GASjoint1, y = GASjoint2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(p1)

p1 = plotData %>%
  ggplot(aes(x = GASindiv1, y = GASindiv2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(p1)

p2 = plotData %>%
  ggplot(aes(x = ECCAjoint1, y = ECCAjoint2, col = genotype, shape = diet)) + geom_point(size = 2.5) + theme(text=element_text(size = 25))
print(p2)
