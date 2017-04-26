clear
load('ZD_Data_5122_minutes.mat')

global u_total tau
test_episodes = u_total;
test_time = 0:5:(5*size(u_total,2)-5);

Gary_results=cell(9,9);
dataG_episodes=cell(9,9);

s=[9,9];
n=81;

q1=zeros(s);
q2=zeros(s);
L2_error = zeros(s);
Linf_error=zeros(s);
AUC_error_sq=zeros(s);
peak_time_error_sq=zeros(s);
peak_height_error_sq=zeros(s);




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

testing = ~logical(eye(s(1)));

for test=1:9
    for train=1:9
        
        fprintf('Computing error measures for training episode %i, test episode %i\n',train,test);
        
        estbractac=rosencode(train,test);
        
        time_vector_mins = estbractac(:,1)*60;
        start_step = find(time_vector_mins==0,1);
        time_vector_tau = time_vector_mins(start_step:5:end);
        u_star = [zeros(1,6), estbractac(start_step:5:end,2)'];
        eplength = size(u_star,2);

        test_u = u_total(test,1:eplength);
        
        
        [peak_est, peaktime_est] = max(u_star);
        [peak_act, peaktime_act] = max(test_u);

        
%         dataG_episodes{train,test}= struct('full_deconvolved_BrAC', u_star,...
%                                             'actual_error', u_star-test_u,...
%                                             'trained_parameters',[eval(sprintf('beta_5122_%i',train)),eval(sprintf('alpha_5122_%i',train))],...
%                                              'L2_error',{sum((u_star-test_u).^2)},...
%                                              'Linf_error',{max(abs(u_star-test_u))},...
%                                              'AUC_error',{sum(u_star)-sum(test_u)},...
%                                              'peak_time_error',{tau*(peaktime_est-peaktime_act)},...
%                                              'peak_height_error',{peak_est-peak_act});
                                         
        q2(train,test)=eval(sprintf('beta_5122_%i',train));
        q1(train,test)=eval(sprintf('alpha_5122_%i',train));
        L2_error(train,test)=sum((u_star-test_u).^2);
        Linf_error(train,test)=max(abs(u_star-test_u));
        AUC_error_sq(train,test)=(sum(u_star)-sum(test_u)).^2;
        peak_time_error_sq(train,test)=(tau*(peaktime_est-peaktime_act)).^2;
        peak_height_error_sq(train,test)=(peak_est-peak_act).^2;
        
    end
    
end


% Pick out test set
testing_array=false(s);
test_ = @(training,test) all(test~=training);

for test = 1:9
    for tr = 1:s(1)
        testing_array(tr,test)=test_(tr,test) & test~=4;
    end
end


ntest = sum(sum(testing_array));
ntrain = 81 - ntest;

        
        

% Now we need the overall mean and stds
overall_q2 = [mean(q2(:)),std(q2(:))];
overall_q1 = [mean(q1(:)),std(q1(:))];

test_stat = @(meas) [mean(reshape(meas(testing_array),ntest,1)),std(reshape(meas(testing_array),ntest,1))];
test_L2 = test_stat(L2_error);
test_Linf = test_stat(Linf_error);
test_AUC_sq = test_stat(AUC_error_sq);
test_peaktime_sq = test_stat(peak_time_error_sq);
test_peakheight_sq = test_stat(peak_height_error_sq);

train_stat = @(meas) [mean(reshape(meas(~testing_array),ntrain,1)),std(reshape(meas(~testing_array),ntrain,1))];
training_L2 = train_stat(L2_error);
training_Linf = train_stat(Linf_error);
training_AUC_sq = train_stat(AUC_error_sq);
training_peaktime_sq = train_stat(peak_time_error_sq);
training_peakheight_sq = train_stat(peak_height_error_sq);

tr_data = [training_AUC_sq,training_peakheight_sq,training_peaktime_sq,training_L2,training_Linf];
te_data = [test_AUC_sq,test_peakheight_sq,test_peaktime_sq,test_L2,test_Linf];

%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     
%     q2_mean(test) = mean(q2(:,test));
%     q2_sd(test) = std(q2(:,test));
%     q1_mean(test) = mean(q1(:,test));
%     q1_sd(test) = std(q1(:,test));
%     
%     % Compute mean and std of estimates for each TEST episode across all
%     % training episodes
%     
%     % For a given test episode, we get mean data across all training
%     % choices besides ''test''.
%     L2_test_mean(test)=mean(L2_error(testing(test,:),test));
%     L2_test_sd(test)=std(L2_error(testing(test,:),test));
%     Linf_test_mean(test)=mean(Linf_error(testing(test,:),test));
%     Linf_test_sd(test)=std(Linf_error(testing(test,:),test));
%     AUC_sq_test_mean(test)=mean(AUC_error_sq(testing(test,:),test));
%     AUC_sq_test_sd(test)=std(AUC_error_sq(testing(test,:),test));
%     peaktime_sq_test_mean(test)=mean(peak_time_error_sq(testing(test,:),test));
%     peaktime_sq_test_sd(test)=std(peak_time_error_sq(testing(test,:),test));
%     peakheight_sq_test_mean(test)=mean(peak_height_error_sq(testing(test,:),test));
%     peakheight_sq_test_sd(test)=std(peak_height_error_sq(testing(test,:),test));
%     
%     L2_training(test) = L2_error(test,test);
%     Linf_training(test) = Linf_error(test,test);
%     AUC_sq_training(test) = AUC_error_sq(test,test);
%     peaktime_sq_training(test) = peak_time_error_sq(test,test);
%     peakheight_sq_training(test) = peak_height_error_sq(test,test);
% 
%     
% end
% 
% 
% 
% % Plot results
% %
% % for test = 1:9
% %     x = dataG_episodes{1,test};
% %     subplot(3,3,test);
% %     plot(x.full_deconvolved_BrAC-x.actual_error,'k--');
% %     hold on
% %     for train = 1:9
% %         y = dataG_episodes{train,test};
% %         plot(y.full_deconvolved_BrAC);
% %     end
% % end
% 
% n=81;
% % Now we need the overall mean and stds
% overall_q2 = [mean(q2_mean),std(reshape(q2,n,1))];
% overall_q1 = [mean(q1_mean),std(reshape(q1,n,1))];
% test_L2 = [mean(L2_test_mean),std(reshape(L2_error,n,1))];
% test_Linf = [mean(Linf_test_mean),std(reshape(Linf_error,n,1))];
% test_AUC_sq = [mean(AUC_sq_test_mean),std(reshape(AUC_error_sq,n,1))];
% test_peaktime_sq = [mean(peaktime_sq_test_mean),std(reshape(peak_time_error_sq,n,1))];
% test_peakheight_sq = [mean(peakheight_sq_test_mean),std(reshape(peak_height_error_sq,n,1))];
% 
% training_L2 = [mean(L2_training),std(L2_training)];
% training_Linf = [mean(Linf_training),std(Linf_training)];
% training_AUC_sq = [mean(AUC_sq_training),std(AUC_sq_training)];
% training_peaktime_sq = [mean(peaktime_sq_training),std(peaktime_sq_training)];
% training_peakheight_sq = [mean(peakheight_sq_training),std(peakheight_sq_training)];
% 
% tr_data = [training_AUC_sq,training_peakheight_sq,training_peaktime_sq,training_L2,training_Linf];
% te_data = [test_AUC_sq,test_peakheight_sq,test_peaktime_sq,test_L2,test_Linf];
% 
% % Plot results
% %
% % for test = 1:9
% %     x = dataG_episodes{1,test};
% %     subplot(3,3,test);
% %     plot(x.full_deconvolved_BrAC-x.actual_error,'k--');
% %     hold on
% %     for train = 1:9
% %         y = dataG_episodes{train,test};
% %         plot(y.full_deconvolved_BrAC);
% %     end
% % end

save('rosen_results_notestep4.mat');
