%% Set seed
rng(2);

%% Set RKHS
rkhs_eigenfile = 'green1_eigen';
data_path = 'none';
% data_path = '030117_234splhr_fixedtraining_testeps15_arrays';

%% Set hyperparameters
P = 10;   %number of eigenfunctions
T = 1;
n = 50;
tau = 1/50;
t = 0:tau:(T-tau);


%% Generate initial theta

% number of parameters
nTheta=2;

% initial guess for number of iterations
K=60;

% initialize empty arrays to hold parameter values
thetas=zeros(nTheta,K);

% Initial parameter guesses
thetas(1,1)= 1;
thetas(2,1)= 1;


% initialize empty arrays to hold amplitudes
a = zeros(P,K);

%% Get eigenvalues and eigenfunctions
[eivs,eifs] = get_kernel_eigenstuff(P,T,rkhs_eigenfile);


%% Generate sample input
% Base signal plus noise
actual_u = pb_7p2_example_u() ;
sampled_u = actual_u(t);

% plot(t,sampled_u);


%% Generate sample output
% fprintf('Generating sample output y')

% y = zeros(1,n);
% y(1)=0;
% for i = 1:(n-1)
%     y(i+1) = L_i([5,9.22],actual_u,i,tau);
% end

y = [0,7.88095366723674e-05,0.00103125740020056,0.00425802631203408,0.0109487106742877,0.0216727342166110,0.0363299402796804,0.0542128037925836,0.0742872217889489,0.0952438390517745,0.115773735485898,0.134695370184833,0.151113977012157,0.164215096458382,0.173782968335687,0.179509850250412,0.181500188617153,0.180170029832357,0.176002392924988,0.169616154122135,0.161547784539635,0.152696662331240,0.143452966474724,0.134612782969953,0.126611577696363,0.119861421105558,0.114578590320305,0.110958828816849,0.108905374421002,0.108205189637280,0.108760246571796,0.109861127152109,0.111364541302021,0.112570195812181,0.112986485686659,0.112348166641714,0.110225955711135,0.106496590351219,0.101132319193791,0.0942716224384557,0.0861827407805335,0.0773306326961579,0.0680696708969857,0.0589579388905687,0.0502928147264237,0.0424323522743800,0.0355219726659830,0.0295923564096101,0.0246246942156294,0.0204888398013944];


%plot(t,y)

%% Generate parameter variance estimates 

alph=0.7;
Vhat_ = Vhat(y,[exp(1),10],tau,P,T,n,eivs,rkhs_eigenfile,data_path);
Vhat_ = diag(diag(Vhat_));
% Cell array of matrices, one for each parameter.

%% MCMC

% These are parameters for the convergence condition.
epsilon = 10^(-5);
nEndingSteps=5;

k=2;
rejected=0;
while ~converged(a,k,nEndingSteps,epsilon) %CONVERGENCE CONDITION NOT IMPLEMENTED
        
    c = unifrnd(0,1,nTheta,1); % Generates a uniform random number for each parameter
    
    thetas(:,k)=mvnrnd(thetas(:,k-1),alph*Vhat_);
    
    acc = acceptance(thetas(:,k),thetas(:,k-1),y,tau,T,P,n,rkhs_eigenfile,data_path);
    
    for i = 1:nTheta
        if c(i) > acc(i)
            rejected = rejected + 1;
            thetas(i,k)=thetas(i,k-1);
        end
    end

    
    EV=EV_aP_given_theta_y(thetas(:,k),y,rkhs_eigenfile,P,T,n,tau,data_path);
    
    mus = EV{1};
    covs = EV{2};
    
    
    for i = 1:P
        a(i,k-1) = normrnd(mus(i),sqrt(covs(i,i)));
    end
    
    k=k+1;
    
end

1-rejected/k
