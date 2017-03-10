% process output of deterministic algorithm into data we need

load('030117_234splhr_fixedtraining_testeps15_arrays.mat');

qs = cell2mat(permute(trained_parameters(1,1,1,1,:),[5,1:4]));

q2s = qs(:,1);
q1s = q1(:,2);

q2_mu = mean(q2s);
q2_std = std(q2s);
q1_mu = mean(q1s);
q1_std = std(q1s);

% Also somehow process the error terms