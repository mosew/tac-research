%% Make tables to compare my data to Gary's
global u_total tau y_total

%% Now compute error for Gary's

% Change directory to his
cd('C:\Users\mose\Dropbox\research\Filter Design MW')

% Set episode to train filter on
train = 7;

% Load optimal parameters for filter training episode
q_1_in=eval(sprintf('alpha_5122_%i',train));
q_2_in=eval(sprintf('beta_5122_%i',train));
q_3_in=0;

n=240;

if train>3
    train = train-1;
end
if train>9
    train = train-1;
end

% BrAC of training episode
BrAC=[(0:(n-1))/(60/tau); u_total(train,:)]';

% TAC of training episode (not pre-processed)
%TAC=[eval(sprintf('t_TAC_5122_%i',train))', eval(sprintf('data_TAC_5122_%i',train))];

% TAC of training episode (pre-processed)
TAC = [(0:(n-1))/(60/tau); y_total(train,:)]';

% Get filter
[q_1_out,q_2_out,q_3_out,r1_r2_h,Est_TAC_q_in,Est_TAC_q_out] = BrAC_Estimator_Filter_Design_2(q_1_in,q_2_in,q_3_in,BrAC,TAC);

G_actual_error=cell(5,1);
G_peak_height_error=zeros(5,1);
G_peak_time_error=zeros(5,1);

nis = zeros(5,1);

% For each possible test episode
for test = 1:6

    % TAC signal never returns to 0 in episode 3
    if test==3
        continue
    end
    
    % Make array index for test episode
    ti=test;
    % Correct array index if we're past episode 3
    if test>3
        ti = test-1;
    end

    n=240;
    
    % Get sampled BrAC test data
    BrAC_test=[(0:(n-1))/(60/tau); u_total(ti,:)]';
    
    % Get sampled TAC test data // maybe use y_total here
    TAC_test=[(0:(n-1))/(60/tau);y_total(ti,:)]';

    % Get Gary's estimated time,BrAC,TAC
    rosen_tBT = BrAC_Est_0_G_1_FD(TAC_test,r1_r2_h);

    % Extract Gary's estimated BrAC
    u_star = rosen_tBT(:,2);
    
    % Get processed (i.e. de-noised and possibly zero-padded) BrAC signal
    test_u = u_total(ti,:);
    
    % Get data about peak time and height
    [peak_est, peaktime_est] = max(u_star);
    [peak_act, peaktime_act] = max(test_u);

    % Set max episode length in timesteps
    n=200;
    
    % Store episode length
    nis(ti) = min(fix(size(rosen_tBT,1)/5),n);
    
    % Get error measures for current test episode
    G_actual_error{ti} = (rosen_tBT(1:5:(5*nis(ti)),2)-u_total(ti,1:nis(ti))');
    G_peak_time_error(ti) = 5*(peaktime_est-peaktime_act);
    G_peak_height_error(ti) = peak_est-peak_act;
    
    % Plot
    figure
    plot(60*rosen_tBT(1:5:(5*nis(ti)),1),rosen_tBT(1:5:(5*nis(ti)),2))
    hold on
    plot(1:5:(5*nis(ti)),u_total(ti,1:nis(ti)))
end

G_L2_MSE = mean(cell2mat(G_actual_error').^2);
G_L2_MSE_sd = std(cell2mat(G_actual_error).^2,0,1);
G_peak_time_MSE = mean(G_peak_time_error.^2);
G_peak_time_MSE_sd = std(G_peak_time_error.^2,0,1);
G_peak_height_MSE = mean(G_peak_height_error.^2);
G_peak_height_MSE_sd = std(G_peak_height_error.^2,0,1);

%% Then me
% Load my own data
load('030117_234splhr_fixedtraining_testeps15_arrays')

n=200;

actual_errors = permute(actual_errors,[2,4,5,1,3]);
% 2 spl/hr
actual_errors = permute(actual_errors(1,:,:),[2,3,1]);


% Arrange computed MSE across the test episodes 1:5 (except 3 and 9 remembeR)
% L2_MSE = permute(L2_MSE,[2,4,1,3]);
% L2_MSE_sd = permute(L2_MSE_sd,[2,4,1,3]);
peak_height_MSE = permute(peak_height_MSE,[2,4,1,3]);
peak_height_MSE_sd = permute(peak_height_MSE_sd,[2,4,1,3]);
peak_time_MSE = permute(peak_time_MSE,[2,4,1,3]);
peak_time_MSE_sd = permute(peak_time_MSE_sd,[2,4,1,3]);


% We just want 2 splines per hour
% L2_MSE = permute(L2_MSE(1,:),[2,1]);
% L2_MSE_sd = permute(L2_MSE_sd(1,:),[2,1]);
peak_height_MSE = permute(peak_height_MSE(1,:),[2,1]);
peak_height_MSE_sd = permute(peak_height_MSE_sd(1,:),[2,1]);
peak_time_MSE = permute(peak_time_MSE(1,:),[2,1]);
peak_time_MSE_sd = permute(peak_time_MSE_sd(1,:),[2,1]);


L2_errors = cell2mat(actual_errors(1,:)').^2;
L2_MSE = mean(L2_errors(:,1:200),2);
L2_MSE_sd = std(L2_errors(:,1:200),0,2);


[G_L2_MSE',L2_MSE]
[G_peak_height_MSE,G_peak_height_MSE_sd;peak_height_MSE,peak_height_MSE_sd]
[G_peak_time_MSE,G_peak_time_MSE_sd;peak_time_MSE,peak_time_MSE_sd]
