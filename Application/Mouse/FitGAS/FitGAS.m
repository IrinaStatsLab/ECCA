clear all
load('Mouse.mat')
%x1 follows normal, with dimensions 40 by 120;
%x2 follows binomial, with dimensions 40 by 21.

distr1 = 'normal';
distr2 = 'binomial';
% Setup 1 -- same as ECCA, r0 = 2, r1 = 1, r2 = 2;
[U0_set1,U1_set1,U2_set1,V0_set1,A1_set1,A2_set1,Mean0_set1,flag_set1] = GAS(X1,X2,2,1,2,distr1,distr2);

% Setup 2: r0 = 1, r1 = 2, r2 = 3;
[U0_set2,U1_set2,U2_set2,V0_set2,A1_set2,A2_set2,Mean0_set2,flag_set2] = GAS(X1,X2,1,2,3,distr1,distr2);

save('Mice_GAS_output.mat','U0_set1', 'U1_set1', 'U2_set1', 'V0_set1', 'A1_set1', 'A2_set1', 'Mean0_set1', 'flag_set1', 'U0_set2', 'U1_set2', 'U2_set2', 'V0_set2', 'A1_set2', 'A2_set2', 'Mean0_set2', 'flag_set2');

saveR('../GAS_result.RData', 'U0_set1', 'U1_set1', 'U2_set1', 'V0_set1', 'A1_set1', 'A2_set1', 'Mean0_set1', 'flag_set1', 'U0_set2', 'U1_set2', 'U2_set2', 'V0_set2', 'A1_set2', 'A2_set2', 'Mean0_set2', 'flag_set2')