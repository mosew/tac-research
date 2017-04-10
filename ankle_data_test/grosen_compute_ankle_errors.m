%% RUN GARY'S DECONVOLUTION ON ANKLE DATA AND COLLECT ERROR MEASURES


train = 0;

% Load optimal parameters for filter training episode
q_1_in=eval(sprintf('alpha_5122_%i',train));
q_2_in=eval(sprintf('beta_5122_%i',train));
q_3_in=0;

% Number of TIMESTEPS per episode
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
G_AUC_error=zeros(5,1);

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
    
    t=(0:(n-1))/(60/tau);
    % Get sampled BrAC test data. First column is time vector, second is data
    BrAC_test=[t; u_total(ti,:)]';
    
    % Get sampled (processed) TAC test data
    TAC_test=[t;y_total(ti,:)]';

    % Get Gary's estimated time,BrAC,TAC
    rosen_tBT = BrAC_Est_0_G_1_FD(TAC_test,r1_r2_h);

    % Extract Gary's estimated BrAC
    u_star = rosen_tBT(:,2);
    
    % Get preprocessed (i.e. de-noised and possibly zero-padded) BrAC signal
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
    G_peak_time_error(ti) = (peaktime_est-peaktime_act*5);
    G_peak_height_error(ti) = peak_est-peak_act;
    G_AUC_error(ti) = sum(rosen_tBT(1:(5*nis(ti)),2))-sum(u_total(ti,1:nis(ti))')*5;
    
    % Plot
    figure
    % Plot Gary's estimated BrAC
    plot(60*rosen_tBT(1:5:(5*nis(ti)),1),rosen_tBT(1:5:(5*nis(ti)),2))
    hold on
