% reprocess data
clear
load('ZD_Data_5122.mat');
data_TAC_5122 = cell(1,11);
t_TAC_5122 = cell(1,11);
for i = 1:11
    data_TAC_5122{i} = eval(sprintf('data_TAC_5122_%d',i));
    t_TAC_5122{i} = eval(sprintf('t_TAC_5122_%d',i))*60;
    t_BrAC_5122{i} = t_BrAC_5122{i}*60;
end
save('ZD_Data_5122_minutes.mat');