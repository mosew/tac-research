%% Set seed
rng(3);

%% Set RKHS
rkhs_eigenfile = 'green1_eigen';
% data_path = 'none';
data_path = 'my_results';

%% Set hyperparameters
test_ep = 2;
N = 32;   %space discretization
P = 6;   %number of eigenfunctions
tau = 5;
burnin = 1;
K = 3;
alph = 3333; % Adjusts the variance of the parameter draw; according to paper should be s.t. acceptance rate is .23


% w = warning('query','last');
% id=w.identifier;
% warning('off',id);

%% Generate initial theta

% number of parameters
nTheta=3;

% initialize empty arrays to hold parameter values
thetas=zeros(nTheta,K);

% Initial parameter guesses
thetas(1,1)= .0044;
thetas(2,1)= 1.23;
thetas(3,1)=1;

% initialize empty arrays to hold amplitudes
a = zeros(P,K-burnin);

%% Compute operators and filter function handle
Kq = build_Kq(theta(1),N);
AN = build_AN(Kq,theta(2),N);
dAN_dq1=build_dAN_dqM(AN);
d2AN_dq1=zeros(size(dAN_dq1));
[ANhat,dANhat_dq1]=build_expm_stuff(AN,tau);
BN=build_BN(theta(2),N);
BNhat=build_BNhat(AN,ANhat,BN);

save('operators.mat','Kq','AN','dAN_dq1','d2AN_dq1','ANhat','dANhat_dqM','BN','BNhat','tau');  
fprintf('Calculated operators...\n');
%% Get eigenvalues and eigenfunctions
[eivs,eifs] = get_kernel_eigenstuff(P,T,rkhs_eigenfile);

%% Get sample input/output
load('my_results.mat','u_total','y_total');
test_u = u_total(test_ep,:);
y = y_total(test_ep,:);
tau = 30;

n=length(test_u);
T = n*tau;
t = tau:tau:T;

% assert(n==length(y));
% fprintf('Loaded sample input and output, plotting...Press enter to continue\n')
% plot(t,test_u);
% hold on
% plot(t,y);
% pause

%% Generate parameter variance estimates 
% I think technically this should be re-evaluated at the current draw of theta,
% but I don't think it should vary all that much.
fprintf('Calculating Vhat, the parameter variance estimate\n');
Vhat_ = Vhat(test_ep,N,y,[.0044,1.23]',tau,P,T,n,eivs,rkhs_eigenfile,data_path);
Vhat_ = (Vhat_ + Vhat_' ) / 2;
save('operators.mat','Vhat_');
fprintf('Calculated Vhat\n');

Vhat_
pause

%% Initialize steps
k=2;
rejected=0;
progress=struct();

%% Loop
while k<K

    k

    % Restricted to be nonnegative
    thetas(:,k) = rmvnrnd(thetas(:,k-1),alph*Vhat_,1,[-1,0;0,-1],[0,0]');
    
    try
        acc = acceptance(thetas(:,k),thetas(:,k-1),y,tau,T,P,n,rkhs_eigenfile,data_path);
    catch
        fprintf('Something went wrong with calculating acceptance ratio for this draw of theta: %d,%d.\n',thetas(1,k),thetas(2,k));
        thetas(:,k)=thetas(:,k-1);
        acc = acceptance(thetas(:,k),thetas(:,k-1),y,tau,T,P,n,rkhs_eigenfile,data_path);
    end
    
    
    c = unifrnd(0,1,1);
    
    if c > acc
        if k>burnin
            rejected = rejected + 1;
        end
        thetas(:,k)=thetas(:,k-1);
    end

    if k<burnin
        k=k+1;
        continue
    end
    
    EV=EV_aP_given_theta_y(thetas(:,k),y,rkhs_eigenfile,P,T,n,tau,data_path);
    
    mus = EV{1};
    covs = EV{2};
    
    a(:,k+1-burnin) = mvnrnd(mus,covs);
    
%     progress.step = k;
%     progress.theta1=[thetas(1,k-1),thetas(1,k)];
%     progress.theta2=[thetas(2,k-1),thetas(2,k)];
%     progress.acc1=acc(1);
%     progress.acc2=acc(2);
%     progress
    
    
    k=k+1;
    
end

% MCMC acceptance rate for theta
1-rejected/(k-burnin)

fks = cell(1,K-burnin);
for i = 1:K-burnin
    fks{i} = f_from_a_eifs(a(:,i)',eifs);
end

fL_fU_fM = confidence_limits(fks);
fL = fL_fU_fM{1};
fU = fL_fU_fM{2};
fM = fL_fU_fM{3};


plot(fM(t),'b')
hold on
plot(sampled_u,'ko')
plot(fL(t),'r--')
plot(fU(t),'r--')

% Scatterplot thetas
figure
scatter(thetas(1,burnin:end),thetas(2,burnin:end),'.')