# https://mixomicsteam.github.io/Bookdown/ for the tutorial for mixOmics
rm(list = ls())

library(mixOmics)
?nutrimouse
data(nutrimouse)

# gene: expressions of 120 genes measured in liver cells, 
# selected (among about 30,000) as potentially relevant in the context of the nutrition study. 
# These expressions come from a nylon macroarray with radioactive labelling;
# lipid: concentrations (in percentages) of 21 hepatic fatty acids measured by gas chromatography

X1 = as.matrix(nutrimouse$gene)
X2 = as.matrix(nutrimouse$lipid/100) # Scale the percentages to be in the decimal format.

# Use trials 100 to modify x2, same as ECCA.
source("Function/exp_family.R")
X2 = modifyBinomial(X2, trials = 100)

# The estimated principal angles are: 34.98249 57.18569 74.06153
# In ECCA, we use r0 = 2, r1 = 1, r2 = 2
# In GAS, we restrict the joint score matrix to be the same. We consider two possible setup:
# Setup 1: same as ECCA, r0 = 2, r1 = 1, r2 = 2;
# Setup 2: r0 = 1, r1 = 2, r2 = 3;
library('R.matlab')
library('matlabr')

# Save data as the mat file
writeMat('Application/Mouse/FitGAS/Mouse.mat', X1 = X1, X2 = X2)

# Run GAS and save in R
# Attention: may need to run directly in Matlab as performance of saveR is system and R version dependent
run_matlab_script('Application/Mouse/FitGAS/FitGAS.m')
