load('test_data.mat')
load('test_data_constantm.mat')
t=test_data;
t=permute(t,[2,3,1,4]);
t0=test_data_constantm;
t0=permute(t0,[2,3,1,4]);

badscales=[];
L2_errors = zeros(7,5,3,11);
Linf_errors = zeros(7,5,3,11);
AUC_errors = zeros(7,5,3,11);
peak_time_errors = zeros(7,5,3,11);
peak_height_errors = zeros(7,5,3,11);
actual_errors = cell(7,5,3,11);
full_deconvolved_BrACs = cell(7,5,3,11);
trained_parameters = cell(7,5,3,11);

for k = 1:3
    for i = 1:11
        for m = 1:7
            for p = 1:5
                if m==1
                    s = t0(m,p,k,i);
                else
                    s=t(m-1,p,k,i);
                end
                s=s{1};
                if s.badscale==1
                    badscales=[badscales;[m,p,k,i]];
                end
                L2_errors(m,p,k,i) = s.L2_error;
                Linf_errors(m,p,k,i) = s.Linf_error;
                AUC_errors(m,p,k,i) = s.AUC_error;
                peak_time_errors(m,p,k,i) = s.peak_time_error;
                peak_height_errors(m,p,k,i) = s.peak_height_error;
                actual_errors{m,p,k,i} = s.actual_error;
                full_deconvolved_BrACs{m,p,k,i} = s.full_deconvolved_BrAC;
                trained_parameters{m,p,k,i} = s.trained_parameters;
                if m==1
                    trained_parameters{m,p,k,i} = trained_parameters{m,p,k,i}(1:2);
                end
            end
        end
    end
end

q = t;

% These are the runs where a "Nearly Singular Matrix" error was thrown
badscales_actual_mp = badscales;
badscales_actual_mp(:,1) = 2.^(badscales(:,1)-1);
badscales_actual_mp(:,2) = 2.^(badscales(:,2)+1);


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
