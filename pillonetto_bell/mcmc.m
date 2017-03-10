%% Set seed
rng(2);

%% Set RKHS
rkhs_eigenfile = 'greenkernel1';

%% Generate initial parameter data

% number of parameters
nParams=2;

% initial guess for number of iterations
K=200;

% initialize empty arrays to hold parameter values
params=zeros(nParams,K);

% Initial parameter guesses
params(1,1)= 1;
params(2,1)= 10;


% initialize empty arrays to hold amplitudes
a = zeros(P,K);



%% Generate sample input
% Base signal
actual_u = pb_7p2_example_u(0:1/49:1) ;
% plot(actual_u)
% hold on

% Plus noise
actual_u = actual_u + normrnd(0,.05*actual_u);
% plot(actual_u,'o')

%% Generate parameter variance estimates 

% Each of these should give arrays of size nParams x 1
alph=get_alpha_vel(); % Ordinary array NOT IMPLEMENTED
Vhat=Vhat(2); % Cell array of matrices, one for each parameter.

%% MCMC

% These are parameters for the convergence condition.
epsilon = .001;
nEndingSteps=3;

k=1;
while ~converged(a,k,nEndingSteps,epsilon) %CONVERGENCE CONDITION NOT IMPLEMENTED
    k=k+1;
        
    c = unifrnd(0,1,nParams,1); % Generates a uniform random number for each parameter
    
    for i=1:nParams
        params(i,k)=normrnd(params(i,k-1),alph(i)*Vhat(i));
        if c(i) > acceptance(params(i,k),params(i,k-1),y)
            params(i,k)=params(i,k-1);
        end
    end
    
    EV=EV_aP_given_theta_y(theta,y,rkhs_eigenfile,P,T);
    
    a(:,k) = normrnd(EV{1},EV{2});
end