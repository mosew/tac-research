clear

load('ZD_Data_5122_minutes.mat');
load('0414_2splhr_8training_alltesteps.mat','b');

%b=permute(b(1,2,1,:,:),[4,5,1,2,3]);

n = numel(b);
s = size(b);

q1=zeros(s);
q2=zeros(s);
L2_errors = zeros(s);
Linf_errors=zeros(s);
AUC_errors=zeros(s);
peak_time_errors=zeros(s);
peak_height_errors=zeros(s);


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

for test = 1:9
    for train = 1:size(b,1)
        r = b{train,test};
        parms = r.trained_parameters;
        q2(train,test) = parms(1);
        q1(train,test) = parms(2);
        L2_errors(train,test) = r.L2_error;
        Linf_errors(train,test)= r.Linf_error;
        AUC_errors(train,test)= r.AUC_error;
        peak_time_errors(train,test) = r.peak_time_error;
        peak_height_errors(train,test)=r.peak_height_error;
    end

    q2_mean(test) = mean(q2(:,test));
    q2_sd(test) = std(q2(:,test));
    q1_mean(test)=mean(q1(:,test));
    q1_sd(test)=std(q1(:,test));

    L2_mean(test)=mean(L2_errors(:,test));
    L2_sd(test)=std(L2_errors(:,test));
    Linf_mean(test)=mean(Linf_errors(:,test));
    Linf_sd(test)=std(Linf_errors(:,test));
    AUC_mean(test)=mean(AUC_errors(:,test));
    AUC_sd(test)=std(AUC_errors(:,test));
    peaktime_mean(test)=mean(peak_time_errors(:,test));
    peaktime_sd(test)=std(peak_time_errors(:,test));
    peakheight_mean(test)=mean(peak_height_errors(:,test));
    peakheight_sd(test)=std(peak_height_errors(:,test));

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
overall_L2 = [mean(L2_mean),std(reshape(L2_errors,n,1))];
overall_Linf = [mean(Linf_mean),std(reshape(Linf_errors,n,1))];
overall_AUC = [mean(AUC_mean),std(reshape(AUC_errors,n,1))];
overall_peaktime = [mean(peaktime_mean),std(reshape(peak_time_errors,n,1))];
overall_peakheight = [mean(peakheight_mean),std(reshape(peak_height_errors,n,1))];

save('my_results_8training_2splhr.mat');
