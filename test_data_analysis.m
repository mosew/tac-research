


% m, p, lambda, paradigm, episode

bsize = size(b);
badscales=[];
L2_errors = zeros(bsize);
Linf_errors = zeros(bsize);
AUC_errors = zeros(bsize);
peak_time_errors = zeros(bsize);
peak_height_errors = zeros(bsize);
actual_errors = cell(bsize);
full_deconvolved_BrACs = cell(bsize);
trained_parameters = cell(bsize);

for m = 1:3
    for p = 1:2
        for lamb = 1:3
            for para = 1:3
                for i = 1:5
                    s=b(m,p,lamb,para,i);
                    s=s{1};
                    
                    if s.badscale==1
                        badscales=[badscales;[m,p,lamb,para,i]];
                    end
                    
                    L2_errors(m,p,lamb,para,i) = s.L2_error;
                    Linf_errors(m,p,lamb,para,i) = s.Linf_error;
                    AUC_errors(m,p,lamb,para,i) = s.AUC_error;
                    peak_time_errors(m,p,lamb,para,i) = s.peak_time_error;
                    peak_height_errors(m,p,lamb,para,i) = s.peak_height_error;
                    actual_errors{m,p,lamb,para,i} = s.actual_error;
                    full_deconvolved_BrACs{m,p,lamb,para,i} = s.full_deconvolved_BrAC;
                    trained_parameters{m,p,lamb,para,i} = s.trained_parameters;
                end
            end
        end
    end
end


logL2_errors = log(L2_errors);
logL2_MSE = mean(logL2_errors,5);

% Calculate the mean and SD of error measures across test episodes.

L2_MSE = mean(L2_errors,5);
Linf_error_means = mean(Linf_errors,5);
AUC_MSE = mean(AUC_errors.^2,5);
peak_time_MSE = mean(peak_time_errors.^2,5);
peak_height_MSE = mean(peak_height_errors.^2,5);

AUC_error_means = mean(AUC_errors,5);
peak_height_error_means = mean(peak_height_errors,5);
peak_time_error_means = mean(peak_time_errors,5);

L2_MSE_sd = std(L2_errors,0,5);
logL2_sd = std(logL2_errors,0,5);
Linf_sd = std(Linf_errors,0,5);
AUC_error_sd = std(AUC_errors,0,5);
peak_time_sd = std(peak_time_errors,0,5);
peak_height_sd = std(peak_height_errors,0,5);

AUC_MSE_sd = std(AUC_errors.^2,0,5);
peak_time_MSE_sd = std(peak_time_errors.^2,0,5);
peak_height_MSE_sd = std(peak_height_errors.^2,0,5);

['Saving arrays']
save('bigtest_regH012_retry_arrays.mat');

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
