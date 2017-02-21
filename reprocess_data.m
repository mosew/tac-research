%% Load original data
load('ZD_Data_5122.mat');

%% Set globals
global u_total y_total
temp_model=0;

global tau
if isempty(tau)
    tau=5;
end


%% Fix some weirdness in translation.
t_BrAC_5122{8} = [0,0.5,t_BrAC_5122{8}(2:end)];
data_BrAC_5122{8} = [0,32,data_BrAC_5122{8}(2:end)']';
% From the spreadsheet, this 22 was read as a 32, but Susan thought it
% should be 22, for some reason.

%% Pad pads the beginning of each episode with 5*pad minutes of zeros
pad = 6;
cutoff = 2;

data_TAC_5122 = cell(1,11);
t_TAC_5122 = cell(1,11);
durations = zeros(1,11);

for i = 1:11
    
    % Are we correcting for temperature?
    if temp_model ==0 
        data_TAC_5122{i} = [zeros(1,pad),eval(sprintf('data_TAC_5122_%d',i))']; %no temp correction
    else
        if temp_model==1
            data_TAC_5122{i} = max([zeros(1,pad),cell2mat(data_corrected_TAC_5122(i))'],0); %linear temperature-corrected
        else
            print(temp_model);
            return
        end
    end
    
    
    t_TAC_5122{i} = [5*(-pad:-1),eval(sprintf('t_TAC_5122_%d',i))*60]+pad*5;
    t_BrAC_5122{i} = [5*(-pad:-1),t_BrAC_5122{i}*60]+pad*5;
    data_BrAC_5122{i} = [zeros(1,pad),data_BrAC_5122{i}'];
    
    
    % Make "durations" vector where the episode is over if pad consecutive
    % readings <1 are made.
    iftacdata = data_TAC_5122{i}<1;
    for j = 1:(length(data_TAC_5122{i})-cutoff)
        if all(iftacdata(j:j+cutoff) == [1,zeros(1,cutoff)])
             durations(i) = ceil((j+cutoff)*5/tau);
             break
        end
    end        
end


%% We want to exclude episodes 3 and 9
% Episode 3 the TAC signal never returns to 0
% Episode 9 Susan takes off her sensor during the drinking episode    
dontuse = [3,9];
toUseIndex = true(1,11);
toUseIndex(dontuse)=false;
t_TAC_5122 = t_TAC_5122(toUseIndex);
t_BrAC_5122 = t_BrAC_5122(toUseIndex);
data_TAC_5122 = data_TAC_5122(toUseIndex);
data_BrAC_5122 = data_BrAC_5122(toUseIndex);

%% Prepare data
[~,u_total,y_total] = prepare_data(t_TAC_5122,t_BrAC_5122,data_TAC_5122,data_BrAC_5122,5,tau);



%% Temp model?
if temp_model ==0
    if tau~=5
        save(sprintf('ZD_Data_5122_minutes_tau%d.mat',tau));
    else
        save('ZD_Data_5122_minutes.mat');
    end
else
    if temp_model ==1
        save('ZD_Data_5122_minutes_lineartemp_abs.mat')
    else
        print(temp_model);
        return
    end
end