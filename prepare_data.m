function [time_out,u_total,y_total] = prepare_data(t_TAC,t_BrAC,data_TAC,data_BrAC,tau_in,tau_out)

% INPUT:
% 
% t_TAC     : cell array consisting of time vectors for all TAC episodes
% t_BrAC    :              ""                               BrAC   ""
% data_TAC  :              ""          data       ""        TAC    ""
% data_BrAC :              ""                               BrAC   ""
% tau_in    : length (in HOURS) of one unit in input time vectors
% tau_out   : desired output timestep (in HOURS)
%
% Training and test data should be grouped together.
% Time indexing must be consistent across inputs!
% Recommended is to go through all episodes and trim them.
% Also, each episode must start at time index 0 (NOT 1)


% OUTPUT
% time_out  : vector of time indices used to sample splines, with timestep tau_out
% u_total   : ""                data                BrAC    ""   , fit to spline and sampled uniformly with timestep tau
% y_total   : ""                                     TAC    ""




% Total number of episodes
m_total = length(t_TAC);

% Check to make sure time index of every episode starts at t=0
for i=1:m_total
    assert(t_TAC{i}(1)==0)
    assert(t_BrAC{i}(1)==0)
end


% Find maximum episode length (in tau_in timesteps)
max_n_in=0;
for i=1:length(t_TAC)
    max_n_in = max(max(length(t_TAC{i}),length(t_BrAC{i})),max_n_in);
end

% Number of timesteps in output vector
n_out = floor(max_n_in * tau_in / tau_out);


time_out = (0:(n_out-1))*tau_out;


u_total=zeros(m_total,n_out);
y_total=zeros(m_total,n_out);

% Generate, evaluate, restrict>=0 splines.
for i = 1:m_total
    y_total(i,:)=max(interp1(t_TAC{i},data_TAC{i},time_out,'spline','extrap'),0);
    u_total(i,:)=max(interp1(t_BrAC{i},data_BrAC{i},time_out,'linear','extrap'),0);
    %u_total(i,:) = max( ppval(spline(t_BrAC{i},data_BrAC{i}),time_out), 0);
    %y_total(i,:) = max( ppval(spline(t_TAC{i},data_TAC{i}),time_out), 0);
 end
 
% % De-noise TAC spline
% 

y_total = max(sgolayfilt(y_total',2,41),0)';

% y_total_ma = y_total;
% b=2;
% global n
% for i=b+1:n-b
%      y_total_ma(:,i)=mean(y_total(:,i-b:(i+b)),2);
% end
% y_total=y_total_ma;

% order=2;
% framelen=11;
% B = sgolay(order,framelen);
% for i = 1:size(y_total,1)
%     y_total(i,:) = conv(y_total(i,:),B((framelen+1)/2,:),'same');
% end

% % De-noise BrAC spline
%u_total = max(sgolayfilt(u_total',2,15),0)';
