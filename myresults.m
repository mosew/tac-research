clear

load('ZD_Data_5122_minutes.mat');
load('0416_26splhr_M01_1training_alltesteps.mat','b');

b=permute(b(2,2,:,:),[4,1,2,3]);

n = numel(b);
s = size(b);

q1=zeros(s);
q2=zeros(s);
L2_error = zeros(s);
Linf_error=zeros(s);
AUC_error_sq=zeros(s);
peak_time_error_sq=zeros(s);
peak_height_error_sq=zeros(s);


testing = ~logical(eye(s(1)));



q2_mean = zeros(1,9);
q2_sd = zeros(1,9);
q1_mean = zeros(1,9);
q1_sd = zeros(1,9);

% To gather data on mean error per test episode.
L2_test_mean = zeros(1,9);
L2_test_sd = zeros(1,9);
Linf_test_mean = zeros(1,9);
Linf_test_sd = zeros(1,9);
AUC_sq_test_mean = zeros(1,9);
AUC_sq_test_sd = zeros(1,9);
peaktime_sq_test_mean = zeros(1,9);
peaktime_sq_test_sd = zeros(1,9);
peakheight_sq_test_mean = zeros(1,9);
peakheight_sq_test_sd = zeros(1,9);

L2_training = zeros(1,9);
Linf_training = zeros(1,9);
AUC_sq_training = zeros(1,9);
peaktime_sq_training = zeros(1,9);
peakheight_sq_training = zeros(1,9);

% To gather training data

for test = 1:9
    
    
    for train = 1:size(b,1)
        
        r = b{train,test};
        parms = r.trained_parameters;
        q2(train,test) = parms(1);
        q1(train,test) = parms(2);
        L2_error(train,test) = r.L2_error;
        Linf_error(train,test)= r.Linf_error;
        AUC_error_sq(train,test)= r.AUC_sq_error;
        peak_time_error_sq(train,test) = r.peak_time_sq_error;
        peak_height_error_sq(train,test)=r.peak_height_sq_error;
    
    end

    
    q2_mean(test) = mean(q2(:,test));
    q2_sd(test) = std(q2(:,test));
    q1_mean(test) = mean(q1(:,test));
    q1_sd(test) = std(q1(:,test));
    
    % Compute mean and std of estimates for each TEST episode across all
    % training episodes
    
    L2_test_mean(test)=mean(L2_error(testing(test,:),test));
    L2_test_sd(test)=std(L2_error(testing(test,:),test));
    Linf_test_mean(test)=mean(Linf_error(testing(test,:),test));
    Linf_test_sd(test)=std(Linf_error(testing(test,:),test));
    AUC_sq_test_mean(test)=mean(AUC_error_sq(testing(test,:),test));
    AUC_sq_test_sd(test)=std(AUC_error_sq(testing(test,:),test));
    peaktime_sq_test_mean(test)=mean(peak_time_error_sq(testing(test,:),test));
    peaktime_sq_test_sd(test)=std(peak_time_error_sq(testing(test,:),test));
    peakheight_sq_test_mean(test)=mean(peak_height_error_sq(testing(test,:),test));
    peakheight_sq_test_sd(test)=std(peak_height_error_sq(testing(test,:),test));
    
    L2_training(test) = L2_error(test,test);
    Linf_training(test) = Linf_error(test,test);
    AUC_sq_training(test) = AUC_error_sq(test,test);
    peaktime_sq_training(test) = peak_time_error_sq(test,test);
    peakheight_sq_training(test) = peak_height_error_sq(test,test);

    
end




% Plot results
%
% for test = 1:9
%     x = dataG_episodes{1,test};
%     subplot(3,3,test);
%     plot(x.full_deconvolved_BrAC-x.actual_error,'k--');
%     hold on
%     for train = 1:9
%         y = dataG_episodes{train,test};
%         plot(y.full_deconvolved_BrAC);
%     end
% end


% Now we need the overall mean and stds
overall_q2 = [mean(q2_mean),std(reshape(q2,n,1))];
overall_q1 = [mean(q1_mean),std(reshape(q1,n,1))];
test_L2 = [mean(L2_test_mean),std(reshape(L2_error,n,1))];
test_Linf = [mean(Linf_test_mean),std(reshape(Linf_error,n,1))];
test_AUC_sq = [mean(AUC_sq_test_mean),std(reshape(AUC_error_sq,n,1))];
test_peaktime_sq = [mean(peaktime_sq_test_mean),std(reshape(peak_time_error_sq,n,1))];
test_peakheight_sq = [mean(peakheight_sq_test_mean),std(reshape(peak_height_error_sq,n,1))];

training_L2 = [mean(L2_training),std(L2_training)];
training_Linf = [mean(Linf_training),std(Linf_training)];
training_AUC_sq = [mean(AUC_sq_training),std(AUC_sq_training)];
training_peaktime_sq = [mean(peaktime_sq_training),std(peaktime_sq_training)];
training_peakheight_sq = [mean(peakheight_sq_training),std(peakheight_sq_training)];


save('my_results_1training_M1_6splhr.mat');
