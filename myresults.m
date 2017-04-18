clear

load('ZD_Data_5122_minutes.mat');
load('0416_2splhr_8training_alltesteps.mat','b');

%b=permute(b(2,2,:,:),[3,4,1,2]);

n = numel(b);
s = size(b);


q1=zeros(s);
q2=zeros(s);
L2_error = zeros(s);
Linf_error=zeros(s);
AUC_error_sq=zeros(s);
peak_time_error_sq=zeros(s);
peak_height_error_sq=zeros(s);


for trainInd = 1:size(b,1)

    for test = 1:9
            
        r = b{trainInd,test};
        parms = r.trained_parameters;
        q2(trainInd,test) = parms(1);
        q1(trainInd,test) = parms(2);
        L2_error(trainInd,test) = r.L2_error;
        Linf_error(trainInd,test)= r.Linf_error;
        AUC_error_sq(trainInd,test)= r.AUC_sq_error;
        peak_time_error_sq(trainInd,test) = r.peak_time_sq_error;
        peak_height_error_sq(trainInd,test)=r.peak_height_sq_error;        
        
    end

end




% Pick out test set
testing_array=false(s);
test_ = @(training,test) all(test~=training);

for test = 1:9
    for tr = 1:s(1)
        testing_array(tr,test)=test_(tr,test);
    end
end

ntest = sum(sum(testing_array));
ntrain = n - ntest;

        
        

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

save('my_results_8training_2splhr.mat');

