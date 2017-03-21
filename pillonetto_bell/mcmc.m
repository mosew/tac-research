%% Set seed
rng(2);

%% Set RKHS
rkhs_eigenfile = 'green1_eigen';
data_path = 'none';
% data_path = '030117_234splhr_fixedtraining_testeps15_arrays';

%% Set hyperparameters
P = 30;   %number of eigenfunctions
T = 1;
n = 50;
tau = 1/49;
t = 0:tau:T;


%% Generate initial theta

% number of parameters
nTheta=2;

% initial guess for number of iterations
K=25;

% initialize empty arrays to hold parameter values
thetas=zeros(nTheta,K);

% Initial parameter guesses
thetas(1,1)= 1;
thetas(2,1)= 10;


% initialize empty arrays to hold amplitudes
a = zeros(P,K);

%% Get eigenvalues and eigenfunctions
[eivs,eifs] = get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);


%% Generate sample input
% Base signal plus noise
actual_u = pb_7p2_example_u() ;
sampled_u = actual_u(t);

plot(t,sampled_u);


%% Generate sample output
y = zeros(1,n);
fprintf('Generating sample output y')
y(1)=0;
for i = 1:(n-1)
    y(i+1) = L_i([5,9.22],actual_u,i,tau);
end

%% Generate parameter variance estimates 

alph=.02;
Vhat_ = Vhat(y,[5,10],tau,P,T,n,eivs,rkhs_eigenfile,data_path);
% Cell array of matrices, one for each parameter.

%% MCMC

% These are parameters for the convergence condition.
epsilon = 10^(-3);
nEndingSteps=3;

k=1;
rejected=0;
while ~converged(a,k,nEndingSteps,epsilon) %CONVERGENCE CONDITION NOT IMPLEMENTED
    k=k+1;
        
    c = unifrnd(0,1,nTheta,1); % Generates a uniform random number for each parameter
    
    
    for i=1:nTheta
        thetas(i,k)=normrnd(thetas(i,k-1),alph*Vhat_(i));
    end
    
    acc = acceptance(thetas(:,k),thetas(:,k-1),y,tau,T,P,n,rkhs_eigenfile,data_path);
    
    for i = 1:nTheta
        if c(i) > acc(i)
            rejected = rejected + 1;
            thetas(i,k)=thetas(i,k-1);
        end
    end

    
    EV=EV_aP_given_theta_y(theta,y,rkhs_eigenfile,P,T,n,tau,data_path);
    
    mus = EV{1};
    covs = EV{2};
    
    for i = 1:nTheta
        a(i,k) = normrnd(mus(i),sqrt(covs(i,i)));
    end
end

rejected/k
