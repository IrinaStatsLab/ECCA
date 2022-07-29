clear all
load('Mouse.mat')
%x1 follows normal, with dimensions 40 by 120;
%x2 follows binomial, with dimensions 40 by 21.

rcand1 = 1:10;
rcand2 = 1:10;
Nfold = 10;
[temp_avgCVscore1,r1,temp_allCVscore1]=Nfold_CV_Single(x1, 'normal', rcand1, Nfold);
[temp_avgCVscore2,r2,temp_allCVscore2]=Nfold_CV_Single(x2, 'binomial', rcand2, Nfold, struct('Niter',500, 'lambda', 1e-2));
saveR('MouseRank.RData', 'r1', 'r2', 'Nfold', 'temp_allCVscore1', 'temp_allCVscore2')