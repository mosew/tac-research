
load('ZD_Data_5122_minutes.mat')

global tau M N P
    tau=5;% timestep (in MINUTES)
    N=32; % state and evolution operator discretization index
    M=1; % Q discretization index
    P=25; % U discretization index.  is the number of splines per hour.

% Set training and test sets
% global training test total
training = 5;
test = 5;
total=5; %used for setting max time window. needs to be a RANGE that includes training and test indices.

[~,u_total,y_total] = prepare_data(t_TAC_5122(total),t_BrAC_5122(total),data_TAC_5122(total),data_BrAC_5122(total),5,tau);

training = training-total(1)+1;
test = test-total(1)+1;

y_total = max(y_total,0);


global training_u training_y test_y test_u m n
    training_y = y_total(training,:);
    [m,n] = size(training_y);
    training_y = [zeros(m,1),y_total(training,:)];
    training_u = [u_total(training,:),zeros(m,1)];
    test_y = [0,y_total(test,:)];
    test_u = [u_total(test,:),0];
    
    % training_y should be an ordinary array with one ROW per training episode, n columns
    assert(all(size(training_y)==[m n+1]));
    % training_u should be an ordinary array with one ROW per training episode, n columns
    assert(all(size(training_u)==[m n+1]));
    % test_y should be an ordinary 1xn array
    assert(all(size(test_y)==[1 n+1]));


global Y        
    Y=[training_y;test_y];

    
global MSpl L R
    % LINEAR SPLINE MATRIX
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    L = diag([1,zeros(1,N)]);
    R = diag([zeros(1,N),1]);
    
SplineHandles = cell(1,M+1);
for k=0:M
    SplineHandles{k+1} = @(x) (x>(k-1)/M).*(x<=k/M).*(M*x-(k-1)) + (x>k/M).*(x<(k+1)/M).*(-M*x+k+1);
end

global SplinesP_linear SplinesP_gaussian
t=0:1/n:1;
SplinesP_linear = zeros(P+1,n+1);
SplinesP_gaussian = zeros(P+1,n+1);

for k=0:P
    linearfun = @(x)(x>(k-1)/P).*(x<=k/P).*(P*x-(k-1)) + (x>k/P).*(x<(k+1)/P).*(-P*x+k+1);
    SplinesP_linear(k+1,:) = linearfun(t);
    gaussianfun = @(x) exp(-(x-k/P).^2/(2*sqrt(P)));
    SplinesP_gaussian(k+1,:) = gaussianfun(t);
end

global si_diag si_offdiag
    % si short for "spline integrals"
    % each, when multiplied on the left by q1M, gives the diag, lower,
    % upper of Kq, respectively.
    si_diag=zeros(M+1,N+1); % at (k,j) gives integral of psi_k over support of psi_j
    si_offdiag =zeros(M+1,N);
    for k=1:(M+1)
        for j=1:N
            si_diag(k,j)=integral(SplineHandles{k},max(j-2,0)/N,min(j,N)/N);
            si_offdiag(k,j)=integral(SplineHandles{k},max(j-1,0)/N,j/N);
        end
        si_diag(k,N+1)=integral(SplineHandles{k},(N-1)/N,1);
    end
    
    
global dAN_dqM
    dAN_dqM = zeros(N+1,N+1,M+2);
    % First page is q2, all zeros.
    for k=1:M+1
        dAN_dqM(:,:,k+1) = -MSpl\(N^2*(diag(-si_offdiag(k,:),-1)+diag(si_diag(k,:))+diag(-si_offdiag(k,:),1)));
    end
    
        
global CNhat
    CNhat = [1,zeros(1,N)];

    
   
% Regularization
    
lambda = 0.007;
lambda2 = 0.001;

global Reg dReg
    Reg = @(qM,u) lambda*sum(u.^2) + lambda2*sum(diff(u).^2);
    dReg = @(qM,u) [zeros(1,M+2), 2*(lambda*u + lambda2*[diff(u),-u(n)])];
