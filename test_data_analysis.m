%load('test_data_tau5_noeps39_lambda0001.mat')
t=test_data_tau5_noeps39_lambda0001;
% m, p, paradigm, episode

badscales=[];
L2_errors = zeros(5,9);
Linf_errors = zeros(5,9);
AUC_errors = zeros(5,9);
peak_time_errors = zeros(5,9);
peak_height_errors = zeros(5,9);
actual_errors = cell(5,9);
full_deconvolved_BrACs = cell(5,9);
trained_parameters = cell(5,9);

for para = 1:1
    for i = 1:9
        for m = 1:1
            for p = 1:5
%                 if m==1
%                     s = t0(p,k,i);
%                 else
                    s=t(m,p,para,i);
%                 end
                s=s{1};
                
%                 if s.badscale==1
%                     badscales=[badscales;[m,p,para,i]];
%                 end
                L2_errors(m,p,para,i) = s.L2_error;
                Linf_errors(m,p,para,i) = s.Linf_error;
                AUC_errors(m,p,para,i) = s.AUC_error;
                peak_time_errors(m,p,para,i) = s.peak_time_error;
                peak_height_errors(m,p,para,i) = s.peak_height_error;
                actual_errors{m,p,para,i} = s.actual_error;
                full_deconvolved_BrACs{m,p,para,i} = s.full_deconvolved_BrAC;
                trained_parameters{m,p,para,i} = s.trained_parameters;
%                 if m==1
%                     trained_parameters{m,p,para,i} = trained_parameters{m,p,para,i}(1:2);
%                 end
            end
        end
    end
end


logL2_errors = log(L2_errors);
logL2_MSE = mean(logL2_errors,4);

% Calculate the mean and SD of error measures across episodes.
L2_MSE = mean(L2_errors,4);
Linf_error_means = mean(Linf_errors,4);
AUC_MSE = mean(AUC_errors.^2,4);
peak_time_MSE = mean(peak_time_errors.^2,4);
peak_height_MSE = mean(peak_height_errors.^2,4);

AUC_error_means = mean(AUC_errors,4);
peak_height_error_means = mean(peak_height_errors,4);
peak_time_error_means = mean(peak_time_errors,4);

L2_MSE_sd = std(L2_errors,0,4);
logL2_sd = std(logL2_errors,0,4);
Linf_sd = std(Linf_errors,0,4);
AUC_error_sd = std(AUC_errors,0,4);
peak_time_sd = std(peak_time_errors,0,4);
peak_height_sd = std(peak_height_errors,0,4);

AUC_MSE_sd = std(AUC_errors.^2,0,4);
peak_time_MSE_sd = std(peak_time_errors.^2,0,4);
peak_height_MSE_sd = std(peak_height_errors.^2,0,4);

% Calculate mean and standard error across M
L2_MSE_M = mean(L2_errors,1);
Linf_error_means_M = mean(Linf_errors,1);
AUC_MSE_M = mean(AUC_errors.^2,1);
peak_time_MSE_M = mean(peak_time_errors.^2,1);
peak_height_MSE_M = mean(peak_height_errors.^2,1);

AUC_error_means_M = mean(AUC_errors,1);
peak_height_error_means_M = mean(peak_height_errors,1);
peak_time_error_means_M = mean(peak_time_errors,1);


L2_MSE_sd_M = std(L2_errors,0,1);
logL2_sd_M = std(logL2_errors,0,1);
Linf_sd_M = std(Linf_errors,0,1);
AUC_error_sd_M = std(AUC_errors,0,1);
peak_time_sd_M = std(peak_time_errors,0,1);
peak_height_sd_M = std(peak_height_errors,0,1);

AUC_MSE_sd_M = std(AUC_errors.^2,0,1);
peak_time_MSE_sd_M = std(peak_time_errors.^2,0,1);
peak_height_MSE_sd_M = std(peak_height_errors.^2,0,1);
