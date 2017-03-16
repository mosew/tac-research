%% Set seed
rng(2);

%% Set RKHS
rkhs_eigenfile = 'green1_eigen';

%% Set hyperparameters
P = 30;   %number of eigenfunctions
T = 1; %minutes

%% Generate initial theta

% number of parameters
nTheta=2;

% initial guess for number of iterations
K=300;

% initialize empty arrays to hold parameter values
thetas=zeros(nTheta,K);

% Initial parameter guesses
thetas(1,1)= 1;
thetas(2,1)= 10;


% initialize empty arrays to hold amplitudes
a = zeros(P,K);



%% Generate sample input
% Base signal
actual_u = pb_7p2_example_u() ;
% plot(actual_u)
% hold on

% Plus noise
actual_u = @(s) feval(actual_u,s) + normrnd(0,.05*feval(actual_u,s));
% plot(actual_u,'o')
tau = 1/49;


%% Generate sample output
y = zeros(1,50);
for i = 1:50
    y(i) = L_i(9.22,actual_u,i,tau);
end

%% Generate parameter variance estimates 

alph=1;
Vhat = Vhat(theta,tau,T,P,n,rkhs_eigenfile,data_path);
% Cell array of matrices, one for each parameter.

%% MCMC

% These are parameters for the convergence condition.
epsilon = .001;
nEndingSteps=3;

k=1;
rejected=0;
while ~converged(a,k,nEndingSteps,epsilon) %CONVERGENCE CONDITION NOT IMPLEMENTED
    k=k+1;
        
    c = unifrnd(0,1,nTheta,1); % Generates a uniform random number for each parameter
    
    for i=1:nTheta
        thetas(i,k)=normrnd(thetas(i,k-1),alph*Vhat(i));
        if c(i) > acceptance(i,thetas(i,k),thetas(i,k-1),y)
            rejected = rejected + 1;
            thetas(i,k)=thetas(i,k-1);
        end
    end
    
    EV=EV_aP_given_theta_y(theta,y,rkhs_eigenfile,P,T);
    
    a(:,k) = normrnd(EV{1},EV{2})
end

rejected/k
