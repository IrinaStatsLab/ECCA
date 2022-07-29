clear all
load('BBdata.mat')
nrep = 100;
p1 = 30;
p2 = 20;
tempX1 = bigX1(:,1:p1);
tempX2 = bigX2(:,1:p2);
distr1 = 'binomial';
distr2 = 'binomial';
[U0_GAS,U1_GAS,U2_GAS,V0_GAS,A1_GAS,A2_GAS,Mu0_GAS] = GAS(tempX1,tempX2,r0,r1,r2,distr1,distr2);
[U0_GAS_rank1,U1_GAS_rank1,U2_GAS_rank1,V0_GAS_rank1,A1_GAS_rank1,A2_GAS_rank1,Mu0_GAS_rank1] = GAS(tempX1,tempX2,r0_rank1,r1_rank1,r2_rank1,distr1,distr2);

for i = 2:nrep
    disp(i)
    X1 = bigX1(:,((i-1)*p1 + 1):(i * p1));
    X2 = bigX2(:,((i-1)*p2 + 1):(i * p2));
    [tempU0_GAS,tempU1_GAS,tempU2_GAS,tempV0_GAS,tempA1_GAS,tempA2_GAS,tempMu0_GAS] = GAS(X1,X2,r0,r1,r2,distr1,distr2);
    [tempU0_GAS_rank1,tempU1_GAS_rank1,tempU2_GAS_rank1,tempV0_GAS_rank1,tempA1_GAS_rank1,tempA2_GAS_rank1,tempMu0_GAS_rank1] = GAS(X1,X2,r0_rank1,r1_rank1,r2_rank1,distr1,distr2);
    U0_GAS = [U0_GAS,tempU0_GAS];
    U1_GAS = [U1_GAS,tempU1_GAS];
    U2_GAS = [U2_GAS,tempU2_GAS];
    V0_GAS = [V0_GAS,tempV0_GAS];
    A1_GAS = [A1_GAS,tempA1_GAS];
    A2_GAS = [A2_GAS,tempA2_GAS];
    Mu0_GAS = [Mu0_GAS,tempMu0_GAS];
    
    U0_GAS_rank1 = [U0_GAS_rank1,tempU0_GAS_rank1];
    U1_GAS_rank1 = [U1_GAS_rank1,tempU1_GAS_rank1];
    U2_GAS_rank1 = [U2_GAS_rank1,tempU2_GAS_rank1];
    V0_GAS_rank1 = [V0_GAS_rank1,tempV0_GAS_rank1];
    A1_GAS_rank1 = [A1_GAS_rank1,tempA1_GAS_rank1];
    A2_GAS_rank1 = [A2_GAS_rank1,tempA2_GAS_rank1];
    Mu0_GAS_rank1 = [Mu0_GAS_rank1,tempMu0_GAS_rank1];
end   

saveR('../GAS_BB_results.RData', 'U0_GAS', 'U1_GAS', 'U2_GAS', 'V0_GAS', 'A1_GAS', 'A2_GAS', 'Mu0_GAS', 'U0_GAS_rank1', 'U1_GAS_rank1', 'U2_GAS_rank1', 'V0_GAS_rank1', 'A1_GAS_rank1', 'A2_GAS_rank1', 'Mu0_GAS_rank1')