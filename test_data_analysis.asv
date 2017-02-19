load('REGTEST_test_data_tau5_noeps39_lambda_reg01_dreg1.mat')
t=REGTEST_test_data_tau5_noeps39_lambda_reg01_dreg1;
% m, p, paradigm, episode
%t=permute(t,[2,1,3]);
% t0=test_data_constantm;
% t0=permute(t0,[2,3,1,4]);

badscales=[];
L2_errors = zeros(5,1,9);
Linf_errors = zeros(5,1,9);
AUC_errors = zeros(5,1,9);
peak_time_errors = zeros(5,1,9);
peak_height_errors = zeros(5,1,9);
actual_errors = cell(5,1,9);
full_deconvolved_BrACs = cell(5,1,9);
trained_parameters = cell(5,1,9);

for para = 1:1
    for i = 1:9
        %for m = 1:5
            for p = 1:5
%                 if m==1
%                     s = t0(p,k,i);
%                 else
                    s=t(p,para,i);
%                 end
                s=s{1};
                
                
                
                if s.badscale==1
                    badscales=[badscales;[p,para,i]];
                end
                L2_errors(p,para,i) = s.L2_error;
                Linf_errors(p,para,i) = s.Linf_error;
                AUC_errors(p,para,i) = s.AUC_error;
                peak_time_errors(p,para,i) = s.peak_time_error;
                peak_height_errors(p,para,i) = s.peak_height_error;
                actual_errors{p,para,i} = s.actual_error;
                full_deconvolved_BrACs{p,para,i} = s.full_deconvolved_BrAC;
                trained_parameters{p,para,i} = s.trained_parameters;
%                 if m==1
%                     trained_parameters{p,para,i} = trained_parameters{p,para,i}(1:2);
%                 end
            end
        %end
    end
end

% % These are the runs where a "Nearly Singular Matrix" error was thrown
% badscales_actual_mp = badscales;
% badscales_actual_mp(:,1) = 2.^(badscales(:,1)-1);
% badscales_actual_mp(:,2) = 2.^(badscales(:,2)+1);


logL2_errors = log(L2_errors);
logL2_MSE = mean(logL2_errors,3);

% Calculate the mean and SD of error measures across episodes.
L2_MSE = mean(L2_errors,3);
Linf_error_means = mean(Linf_errors,3);
AUC_MSE = mean(AUC_errors.^2,3);
peak_time_MSE = mean(peak_time_errors.^2,3);
peak_height_MSE = mean(peak_height_errors.^2,3);

AUC_error_means = mean(AUC_errors,3);
peak_height_error_means = mean(peak_height_errors,3);
peak_time_error_means = mean(peak_time_errors,3);

L2_MSE_sd = std(L2_errors,0,3);
logL2_sd = std(logL2_errors,0,3);
Linf_sd = std(Linf_errors,0,3);
AUC_error_sd = std(AUC_errors,0,3);
peak_time_sd = std(peak_time_errors,0,3);
peak_height_sd = std(peak_height_errors,0,3);

AUC_MSE_sd = std(AUC_errors.^2,0,3);
peak_time_MSE_sd = std(peak_time_errors.^2,0,3);
peak_height_MSE_sd = std(peak_height_errors.^2,0,3);

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
