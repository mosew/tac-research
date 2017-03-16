% Collect statistics for how well my stuff does vs. Gary's
%
% I think I've already collected data on how well Gary's does 
%
% We'll use 
clear
load('030117_2splhr_fixedtraining_testeps15_arrays.mat')

% We have to get L2 error, peak height error, peak time error, and AUC error for varying #training episodes

% Index is number of training episodes, starting w/episode 6
a = @(m) permute(m(1,1,1,:),[4,1,2,3]);


stats = [a(peak_height_MSE),a(peak_time_MSE),a(AUC_MSE),a(L2_MSE)]