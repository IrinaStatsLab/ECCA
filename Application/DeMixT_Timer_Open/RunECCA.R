### DO NOT RUN - for illustration purposes only
### To create matched DeMixT/TIMER data - please directly download RNAseq data from https://portal.gdc.cancer.gov/
### and apply DeMixT and TIMER on that data
load("prostate_timer_pit3.RData")
require(pracma)
source('Function/exponentialCCA.R')
source('Function/optU.R')
source('Function/SOC.R')
source('Function/newton_method.R')
source('Function/exp_family.R')
dim(its)
dim(timer)
summary(its)
summary(timer)
 

# result3 is the output used in paper based on the rank selection procedure described.
result3 = ecca_v2(x1 = its, x2 = timer, r0 = 1, r1 = 1, r2 = 2,
                  family1 = "binomial", family2 = "binomial",
                  main_effect = FALSE, gamma = 1000, trials = 1)


save(result3, file = "Application/ItsTimer/ecca_result_without_main_effect.RData")
