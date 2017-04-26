clear
load('ZD_Data_5122_minutes.mat');
names = {'26splhr_1training_M01_alltesteps','2splhr_2training_alltesteps','2splhr_4training_alltesteps','2splhr_8training_alltesteps'};

tr_data = zeros(7,10);
te_data = zeros(7,10);
q = zeros(7,4);
outliers = [];
colinds=[];

iii=0;
for w = 1:(length(names)+3)
    if w<=4
        file = names{1};
    else
        file = names{w-3};
    end
    
    load(sprintf('%s.mat',file),'b');
    
    if w <= 4
        z = floor((w+1)/2);
        y = rem(w+1,2)+1;
        b=permute(b(z,y,:,:),[3,4,1,2]);
    end

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
    training_array=false(s);
    
    % Returns true if the episode "test" was not in the training set
    test_ = @(training,test) all(test~=training) & test~=4;
    train_ = @(training,test) any(test==training);

    for test = 1:9
        for tr = 1:s(1)
            testing_array(tr,test)=test_(tr,test);
            training_array(tr,test)=train_(tr,test);
        end
    end

    ntest = sum(sum(testing_array));
    ntrain = sum(sum(training_array));




    % Now we need the overall mean and stds
    overall_q2 = [mean(q2(:)),std(q2(:))];
    overall_q1 = [mean(q1(:)),std(q1(:))];

    test_stat = @(meas) [mean(reshape(meas(testing_array),ntest,1)),std(reshape(meas(testing_array),ntest,1))];
    test_L2 = test_stat(L2_error);
    test_Linf = test_stat(Linf_error);
    test_AUC_sq = test_stat(AUC_error_sq);
    test_peaktime_sq = test_stat(peak_time_error_sq);
    test_peakheight_sq = test_stat(peak_height_error_sq);

    train_stat = @(meas) [mean(reshape(meas(training_array),ntrain,1)),std(reshape(meas(training_array),ntrain,1))];
    training_L2 = train_stat(L2_error);
    training_Linf = train_stat(Linf_error);
    training_AUC_sq = train_stat(AUC_error_sq);
    training_peaktime_sq = train_stat(peak_time_error_sq);
    training_peakheight_sq = train_stat(peak_height_error_sq);

    tr_data(w,:) = [training_AUC_sq,training_peakheight_sq,training_peaktime_sq,training_L2,training_Linf];
    te_data(w,:) = [test_AUC_sq,test_peakheight_sq,test_peaktime_sq,test_L2,test_Linf];
    q(w,:)=[overall_q1,overall_q2];
    
% %   Distributions of the parameter estimates (M=1 omitted)
%     subplot(2,7,w)
%     histogram(q1);
%     subplot(2,7,w+7)
%     histogram(q2);
%     
%   Find outliers
    measures = {'q1','q2','L2_error','Linf_error','AUC_error_sq','peak_time_error_sq','peak_height_error_sq'};
    for j = 1:length(measures)
        vv=eval(measures{j});
        col = reshape(vv,numel(vv),1);
        colinds=[];
        for i = 1:numel(vv)
            if mahal(col(i),col)>20
                colinds=[colinds;i];
                [z,x]=ind2sub(size(vv),i);
                outliers=[outliers;w,z,x];
%                 measures{j};
%                 mahal(col(i),col)
%                 vv(z,x)
            end
        end
    end
    
%     What went wrong with test episode 4?
    for j = 1:length(colinds)
        plot(b{colinds(j)}.full_deconvolved_BrAC,'color',[.8,.8,.8]);
        hold on
    end
    plot(u_total(4,:),'ko');
    xlim([0,276])

%    %   Looks like it just frequently overshot the uptake; could have been caused by a single bad measurement.

end

save('my_results.mat');