clear
load('ZD_Data_5122_minutes.mat')

global u_total tau
test_episodes = u_total;
test_time = 0:5:(5*size(u_total,2)-5);

Gary_results=cell(9,9);
dataG_episodes=cell(9,9);

q1=zeros(9,9);
q2=zeros(9,9);
L2_errors = zeros(9,9);
Linf_errors=zeros(9,9);
AUC_sq_errors=zeros(9,9);
peak_time_sq_errors=zeros(9,9);
peak_height_sq_errors=zeros(9,9);


q2_mean = zeros(1,9);
q2_sd = zeros(1,9);
q1_mean = zeros(1,9);
q1_sd = zeros(1,9);
L2_mean = zeros(1,9);
L2_sd = zeros(1,9);
Linf_mean = zeros(1,9);
Linf_sd = zeros(1,9);
AUC_mean = zeros(1,9);
AUC_sd = zeros(1,9);
peaktime_mean = zeros(1,9);
peaktime_sd = zeros(1,9);
peakheight_mean = zeros(1,9);
peakheight_sd = zeros(1,9);




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

        
        dataG_episodes{train,test}= struct('full_deconvolved_BrAC', u_star,...
                                            'actual_error', u_star-test_u,...
                                            'trained_parameters',[eval(sprintf('beta_5122_%i',train)),eval(sprintf('alpha_5122_%i',train))],...
                                             'L2_error',{sum((u_star-test_u).^2)},...
                                             'Linf_error',{max(abs(u_star-test_u))},...
                                             'AUC_error',{sum(u_star)-sum(test_u)},...
                                             'peak_time_error',{tau*(peaktime_est-peaktime_act)},...
                                             'peak_height_error',{peak_est-peak_act});
                                         
        q2(train,test)=eval(sprintf('beta_5122_%i',train));
        q1(train,test)=eval(sprintf('alpha_5122_%i',train));
        L2_errors(train,test)=sum((u_star-test_u).^2);
        Linf_errors(train,test)=max(abs(u_star-test_u));
        AUC_sq_errors(train,test)=(sum(u_star)-sum(test_u)).^2;
        peak_time_sq_errors(train,test)=(tau*(peaktime_est-peaktime_act)).^2;
        peak_height_sq_errors(train,test)=(peak_est-peak_act).^2;
        
    end
    
    q2_mean(test) = mean(q2(:,test));
    q2_sd(test) = std(q2(:,test));
    q1_mean(test)=mean(q1(:,test));
    q1_sd(test)=std(q1(:,test));
    
    L2_mean(test)=mean(L2_errors(:,test));
    L2_sd(test)=std(L2_errors(:,test));
    Linf_mean(test)=mean(Linf_errors(:,test));
    Linf_sd(test)=std(Linf_errors(:,test));
    AUC_mean(test)=mean(AUC_sq_errors(:,test));
    AUC_sd(test)=std(AUC_sq_errors(:,test));
    peaktime_mean(test)=mean(peak_time_sq_errors(:,test));
    peaktime_sd(test)=std(peak_time_sq_errors(:,test));
    peakheight_mean(test)=mean(peak_height_sq_errors(:,test));
    peakheight_sd(test)=std(peak_height_sq_errors(:,test));
    
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
overall_q2 = [mean(q2_mean),std(reshape(q2,81,1))];
overall_q1 = [mean(q1_mean),std(reshape(q1,81,1))];
overall_L2 = [mean(L2_mean),std(reshape(L2_errors,81,1))];
overall_Linf = [mean(Linf_mean),std(reshape(Linf_errors,81,1))];
overall_AUC = [mean(AUC_mean),std(reshape(AUC_sq_errors,81,1))];
overall_peaktime = [mean(peaktime_mean),std(reshape(peak_time_sq_errors,81,1))];
overall_peakheight = [mean(peakheight_mean),std(reshape(peak_height_sq_errors,81,1))];

save('rosen_results.mat');
