%load('REGTEST_tau5_noeps39_lambda_10^-0.823909.mat')
t=b;
% m, p, paradigm, episode
global tau
badscales=[];
L2_errors = zeros(4,9);
Linf_errors = zeros(4,9);
AUC_errors = zeros(4,9);
peak_time_errors = zeros(4,9);
peak_height_errors = zeros(4,9);
actual_errors = cell(4,9);
full_deconvolved_BrACs = cell(4,9);
trained_parameters = cell(4,9);

for para = 1:1
    for i = 1:9
        for m = 1:1
            for p = 1:4
                s=b(p,i);
                s=s{1};
%                 if s.badscale==1
%                     badscales=[badscales;[p,i]];
%                 end
                L2_errors(p,i) = s.L2_error;
                Linf_errors(p,i) = s.Linf_error;
                AUC_errors(p,i) = s.AUC_error;
                peak_time_errors(p,i) = s.peak_time_error;
                peak_height_errors(p,i) = s.peak_height_error;
                actual_errors{p,i} = s.actual_error;
                full_deconvolved_BrACs{p,i} = s.full_deconvolved_BrAC;
                trained_parameters{p,i} = s.trained_parameters;
            end
        end
    end
end


logL2_errors = log(L2_errors);
logL2_MSE = mean(logL2_errors,2);

% Calculate the mean and SD of error measures across episodes.
L2_MSE = mean(L2_errors,2);
Linf_error_means = mean(Linf_errors,2);
AUC_MSE = mean(AUC_errors.^2,2);
peak_time_MSE = mean(peak_time_errors.^2,2);
peak_height_MSE = mean(peak_height_errors.^2,2);

AUC_error_means = mean(AUC_errors,2);
peak_height_error_means = mean(peak_height_errors,2);
peak_time_error_means = mean(peak_time_errors,2);

L2_MSE_sd = std(L2_errors,0,2);
logL2_sd = std(logL2_errors,0,2);
Linf_sd = std(Linf_errors,0,2);
AUC_error_sd = std(AUC_errors,0,2);
peak_time_sd = std(peak_time_errors,0,2);
peak_height_sd = std(peak_height_errors,0,2);

AUC_MSE_sd = std(AUC_errors.^2,0,2);
peak_time_MSE_sd = std(peak_time_errors.^2,0,2);
peak_height_MSE_sd = std(peak_height_errors.^2,0,2);

['Saving arrays']
save(sprintf('REGTEST_data_arrays_tau5_noeps39_lambda_%1.4g.mat',lambda));

% % Calculate mean and standard error across M
% L2_MSE_M = mean(L2_errors,1);
% Linf_error_means_M = mean(Linf_errors,1);
% AUC_MSE_M = mean(AUC_errors.^2,1);
% peak_time_MSE_M = mean(peak_time_errors.^2,1);
% peak_height_MSE_M = mean(peak_height_errors.^2,1);
% 
% AUC_error_means_M = mean(AUC_errors,1);
% peak_height_error_means_M = mean(peak_height_errors,1);
% peak_time_error_means_M = mean(peak_time_errors,1);
% 
% 
% L2_MSE_sd_M = std(L2_errors,0,1);
% logL2_sd_M = std(logL2_errors,0,1);
% Linf_sd_M = std(Linf_errors,0,1);
% AUC_error_sd_M = std(AUC_errors,0,1);
% peak_time_sd_M = std(peak_time_errors,0,1);
% peak_height_sd_M = std(peak_height_errors,0,1);
% 
% AUC_MSE_sd_M = std(AUC_errors.^2,0,1);
% peak_time_MSE_sd_M = std(peak_time_errors.^2,0,1);
% peak_height_MSE_sd_M = std(peak_height_errors.^2,0,1);
