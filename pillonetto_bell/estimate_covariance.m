load('smalltest_24splhr_alleps_regH02_arrays.mat','actual_errors')
% M=0, 2 splines per hour, lambda=0.1, test paradigm 3 (only trained on 4 episodes)
a = cell2mat(actual_errors(1,1,1,1,:));
a = permute(a,[2,5,1,3,4]);
%a = a(1:140,:);
V = cov(a');