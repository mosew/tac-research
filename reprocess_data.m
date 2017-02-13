% process data

load('ZD_Data_5122.mat');
global tau
tau = 5;
pad = 6;

data_TAC_5122 = cell(1,11);
t_TAC_5122 = cell(1,11);
durations=zeros(1,11);

% pad episodes with pad zeros.
for i = 1:11
    data_TAC_5122{i} = [zeros(1,pad),eval(sprintf('data_TAC_5122_%d',i))'];
    t_TAC_5122{i} = [tau*(-pad:-1),eval(sprintf('t_TAC_5122_%d',i))*60]+pad*tau;
    t_BrAC_5122{i} = [tau*(-pad:-1),t_BrAC_5122{i}*60]+pad*tau;
    data_BrAC_5122{i} = [zeros(1,pad),data_BrAC_5122{i}'];
    
    iftacdata = data_TAC_5122{i}<1;
    for j = 1:(length(data_TAC_5122{i})-pad)
        if iftacdata(j:j+pad) == [1,zeros(1,pad)]
            durations(i) = j+pad;
        end
    end        
end
[~,u_total,y_total] = prepare_data(t_TAC_5122,t_BrAC_5122,data_TAC_5122,data_BrAC_5122,5,tau);

save('ZD_Data_5122_minutes.mat');