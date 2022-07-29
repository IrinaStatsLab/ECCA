rm(list = ls())
load("Simulation/GG_Simulation_SNR5/Data.RData")
library('R.matlab')
library('matlabr')
nrep = 100
r0 = 3
r1 = 4
r2 = 3
r0_rank1 = 1
r1_rank1 = 6
r2_rank1 = 5
bigX1 = X1_true[[1]]
bigX2 = X2_true[[1]]
for (i in 2:nrep){
  bigX1 = cbind(bigX1, X1_true[[i]])
  bigX2 = cbind(bigX2, X2_true[[i]])
}
writeMat('Simulation/GG_Simulation_SNR5/GAS/GGdata.mat', bigX1 = bigX1, bigX2 = bigX2, 
         r0 = r0, r1 = r1, r2 = r2, r0_rank1 = r0_rank1, r1_rank1 = r1_rank1, r2_rank1 = r2_rank1)
run_matlab_script('Simulation/GG_Simulation_SNR5/GAS/Run_GAS.m')
