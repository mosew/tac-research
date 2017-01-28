clear
load('5122.mat')

global tau M N
    tau=5; % timestep (in MINUTES)
    N=32; % state and evolution operator discretization index
    M=1; % Q discretization index

% Set training and test sets
% global training test total
 training = 1:11;
 test = 1;
 total=1:11; %used for setting max time window. this needs to be a RANGE that includes the training and test indices.

[t,u_total,y_total] = prepare_data(t_TAC_5122(total),t_BrAC_5122(total),data_TAC_5122(total),data_BrAC_5122(total),5,tau);

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



global Y

    % training_y should be an ordinary array with one ROW per training episode, n columns
    assert(all(size(training_y)==[m n+1]));
    % training_u should be an ordinary array with one ROW per training episode, n columns
    assert(all(size(training_u)==[m n+1]));
    % test_y should be an ordinary 1xn array
    assert(all(size(test_y)==[1 n+1]));
        
    Y=[training_y;test_y];


    
global MSpl L R
    % LINEAR SPLINE MATRIX
    MSpl = 1/(6*N) * gallery('tridiag',ones(1,N),[2,4*ones(1,N-1),2],ones(1,N)); % [<psi^N_i,psi^N_j>]_{ij}
    L = diag([1,zeros(1,N)]);
    R = diag([zeros(1,N),1]);
    

global si_diag si_offdiag
    % These are (M+1)xN (resp. (M+1)x(N+1)) matrices.
    % SplineIntegrals0(k,j) gives int_{(j-1)/N} ^{j/N} psi^M_k(x) dx
    % SplineIntegrals1(k,j) gives int_{(j-2)/N} ^{(j)/N} psi^M_k(x) dx
    % with integral limits constrained to [0,1].

    SplineHandles = cell(1,M+1);
    for k=0:M
        SplineHandles{k+1} = @(x) (x>(k-1)/M).*(x<=k/M).*(M*x-(k-1)) + (x>k/M).*(x<(k+1)/M).*(-M*x+k+1);
    end

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
    % outputs for single timestep
    CNhat = [1,zeros(1,N)];

    
lambda = 0.07;
lambda2 = 0.16;

global Reg dReg
    Reg = @(qM,u) lambda*sum(u.^2) + lambda2*sum(diff(u).^2);
    dReg = @(qM,u) [zeros(1,M+2),2*(lambda*u + lambda2*[diff(u),-u(n)])];

