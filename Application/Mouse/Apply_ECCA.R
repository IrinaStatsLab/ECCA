# This file pre-process the mouse data and runs ECCA with given ranks.

# https://mixomicsteam.github.io/Bookdown/ for the tutorial for mixOmics
rm(list = ls())
source('Function/exponentialCCA.R')
source('Function/optU.R')
source('Function/SOC.R')
source('Function/newton_method.R')
source('Function/exp_family.R')
source('Function/AnalyticalUpdate.R')

library(mixOmics)
?nutrimouse
data(nutrimouse)

# gene: expressions of 120 genes measured in liver cells, 
# selected (among about 30,000) as potentially relevant in the context of the nutrition study. 
# These expressions come from a nylon macroarray with radioactive labelling;
# lipid: concentrations (in percentages) of 21 hepatic fatty acids measured by gas chromatography

x1 = nutrimouse$gene
x2 = nutrimouse$lipid/100 # Scale the percentages to be in the decimal format.
# Probably it's reasonable to assume that the number of trials is 100.
summary(x1)
summary(x2)
count = 0
n = dim(x2)[1]
p2 = dim(x2)[2]

# Get proportions of 0s in x2:
for (i in 1:n){
  for (j in 1:p2){
    if (x2[i,j] == 0){
      count = count + 1
    }
  }
}
print(count/(n * p2))


# r0 = 2 (joint rank), r1 = 1 (individual rank 1), r2 = 2 (individual rank 2) based on CV and principal angles.
r0 = 2
r1 = 1
r2 = 2

family1 = 'gaussian'
family2 = 'binomial'

# For zerosm apply transformation:
# if x = 0 then set it to be 3/8/(100+3/4)
x2_mod100 = modifyBinomial(x2, trials = 100)

# Apply ECCA
eccaMat100 = ecca_v2(x1, x2_mod100, r0, r1, r2, family1 = family1, family2 = family2, trials = 100)
print(eccaMat100$constraints)

save(eccaMat100, file = "Application/Mouse/ECCA_result.RData")