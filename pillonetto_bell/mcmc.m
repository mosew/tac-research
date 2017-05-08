%% Set seed
rng(3);

%% Set RKHS
rkhs_eigenfile = 'green1_eigen';
data_path = 'none';
% data_path = '030117_234splhr_fixedtraining_testeps15_arrays';

%% Set hyperparameters
P = 6;   %number of eigenfunctions
T = 1;
n = 50;
tau = 1/50;
t = tau:tau:T;
burnin = 200;
K = 2000;
alph = 49; % Adjusts the variance of the parameter estimate


% w = warning('query','last');
% id=w.identifier;
% warning('off',id);

%% Generate initial theta

% number of parameters
nTheta=2;

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
% 
% for i = 1:n
%     y(i) = L_i([5,9.22],actual_u,i,tau);
% end

y = [7.88231656865881e-05,0.00103141481842703,0.00425993020130455,0.0109414572812767,0.0216727450029239,0.0363139318337980,0.0542174735415880,0.0742698619534772,0.0952090587956584,0.115750968434917,0.134715310805808,0.151123127706880,0.164276892247834,0.173775284340800,0.179504813304885,0.181458101525460,0.180172458519080,0.175994829843011,0.169610986983633,0.161628554675008,0.152658082895845,0.143436721543320,0.134590958131702,0.126626967590220,0.119876351996740,0.114597745090360,0.110989979388056,0.108897237248321,0.108265218350298,0.108729701660452,0.109909197350094,0.111354787636132,0.112544153090380,0.113050227885335,0.112295064501640,0.110246606493282,0.106490657257262,0.101095550495306,0.0942956034779595,0.0862019988127680,0.0773258488423665,0.0680937375500919,0.0589552679357708,0.0502834884490228,0.0424141229228888,0.0355148857563672,0.0296078491160773,0.0246289666773745,0.0204823531847323,0.0170388615127910];
% y = [0,7.88095366723674e-05,0.00103125740020056,0.00425802631203408,0.0109487106742877,0.0216727342166110,0.0363299402796804,0.0542128037925836,0.0742872217889489,0.0952438390517745,0.115773735485898,0.134695370184833,0.151113977012157,0.164215096458382,0.173782968335687,0.179509850250412,0.181500188617153,0.180170029832357,0.176002392924988,0.169616154122135,0.161547784539635,0.152696662331240,0.143452966474724,0.134612782969953,0.126611577696363,0.119861421105558,0.114578590320305,0.110958828816849,0.108905374421002,0.108205189637280,0.108760246571796,0.109861127152109,0.111364541302021,0.112570195812181,0.112986485686659,0.112348166641714,0.110225955711135,0.106496590351219,0.101132319193791,0.0942716224384557,0.0861827407805335,0.0773306326961579,0.0680696708969857,0.0589579388905687,0.0502928147264237,0.0424323522743800,0.0355219726659830,0.0295923564096101,0.0246246942156294,0.0204888398013944];


%plot(t,y)

%% Generate parameter variance estimates 

Vhat_ = Vhat(y,[exp(1),10],tau,P,T,n,eivs,rkhs_eigenfile,data_path);
Vhat_ = (Vhat_ + Vhat_' ) / 2;

%% MCMC

% Set parameters for the convergence condition.
epsilon = 10^(-4);
nEndingSteps=5;

% Initialize steps
k=2;
rejected=zeros(nTheta,1);
progress=struct();


while ~converged(a,k,K,nEndingSteps,burnin,epsilon) %PAPER'S CONVERGENCE CONDITION NOT IMPLEMENTED
    if ~rem(k,100)
        k
    end

    c = unifrnd(0,1,nTheta,1); % Generates a uniform random number for each parameter
%     Vhat_ = Vhat(y,thetas(:,k-1),tau,P,T,n,eivs,rkhs_eigenfile,data_path);
%     Vhat_ = (Vhat_ + Vhat_' ) / 2;
% 
%     if ~all(eigs(Vhat_)>0)
%         fprintf('Vhat %i',k)
%         eigs(Vhat_)
%         pause
%     end

    thetas(:,k) = mvnrnd(thetas(:,k-1),alph*Vhat_);
        
    try 
        acc = acceptance(thetas(:,k),thetas(:,k-1),y,tau,T,P,n,rkhs_eigenfile,data_path);
    catch
        thetas(:,k) = mvnrnd(thetas(:,k-1),alph*Vhat_);
        acc = acceptance(thetas(:,k),thetas(:,k-1),y,tau,T,P,n,rkhs_eigenfile,data_path);
    end
    
    
    
    for i = 1:nTheta
        if c(i) > acc(i) || thetas(i,k)<=0
%             [i,c(i),acc(i),thetas(i,k),thetas(i,k-1)]
%             pause
            if k>burnin
                rejected(i) = rejected(i) + 1;
%                 fprintf('Rejected theta(%i) at %d \n',i,c(i));
            end
            thetas(i,k)=thetas(i,k-1);
        else
            Vhat_ = Vhat(y,thetas(:,k),tau,P,T,n,eivs,rkhs_eigenfile,data_path);
            Vhat_ = (Vhat_ + Vhat_' ) / 2;
        end
    end

    if k<burnin
%         fprintf('Burnin step %i\n',k);
%         fprintf('theta(1)=%d, theta(2)=%d\n',thetas(1,k),thetas(2,k));
        k=k+1;
        continue
    end
    
    EV=EV_aP_given_theta_y(thetas(:,k-burnin+1),y,rkhs_eigenfile,P,T,n,tau,data_path);
    
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
ones(nTheta,1)-rejected/(k-burnin)

% Plot sample mean input function
f = f_from_a_eifs(mean(a(:,1:k),2)',eifs);
plot(f(t))
hold on
plot(sampled_u)

% Scatterplot thetas
figure
scatter(thetas(1,burnin:end),thetas(2,burnin:end))