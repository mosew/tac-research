% Relies on globals file
% Uses Neumann boundary conditions (different from the paper!!)

function [JN,dJN] = JN_and_dJN_uspline(qu)

% INPUT:
% qu(1) is q2, the gain coefficient
% qu(2:M+2) is qM, the linear spline coefficients of the diffusivity parameter
% qu(M+3:M+4+P) is u, the linear spline coefficients of the input

% OUTPUT:
% JN is (scalar) value of cost function JN evaluated at qu and
% dJN, the gradient vector of JN, an M+2+P - vector

%%%%%

% We have to make a bunch of variables global so that this function, which
% we pass to the optimization toolbox, doesn't take any of the below to be variables to be optimized;
% we could define them here but then they'll be created every
% time J is called.

global N M P m n
global training_u test_y Y Reg dReg
global CNhat SplinesP_linear SplinesP_gaussian

% Process inputs
q2 = qu(1);
%q3 = qM_and_u(2);
% qM is a vector of linear spline coefficients.
qM = qu(2:M+2);
c = qu((M+3):end);
% c is a vector of spline coefficients
assert(length(c)==P+1);
test_u_sampled=sample_u_spline(c,SplinesP_linear);
total_u = [training_u;test_u_sampled];


Kq = build_Kq(qM);
AN = build_AN(Kq,q2);
BN = build_BN(q2);
[ANhat,dANhat_dqM] = build_expm_stuff(AN);
BNhat = build_BNhat(AN,ANhat,BN);
dBNhat_dqM = build_dBNhat_dqM(AN,ANhat,dANhat_dqM,BN);
Phi = build_Phi(ANhat,BNhat,total_u);




% INITIALIZE ETA
eta = zeros(N+1,n+1,m+1);
% goes from t=0 to t=tau*n. At t=0, everything is 0.
for i = 1:m+1
    % episodes in Y are assumed to start at 0. The length of Y therefore is
    % n+1, to accomodate for the output at time 0.
    eta(:,n+1,i) = 2*(CNhat*Phi(:,n+1,i)-Y(i,n+1))*CNhat';
end


% COMPUTE GRADIENT CONTRIBUTIONS
dJN = zeros(1,M+3+P); % one for each component of (qM,u)
JN=0;

for j=n:-1:1
    
    for i=1:m+1
        % Running everything backwards, we begin at timestep n+1.
        eta(:,j,i) = ANhat' * eta(:,j+1,i) + 2*(CNhat*Phi(:,j,i)-Y(i,j))*CNhat';
        for k=1:M+2
            dJN(k) = dJN(k) + eta(:,j+1,i)'*(dANhat_dqM(:,:,k)*Phi(:,j,i) + dBNhat_dqM(:,:,k)*total_u(i,j)); % + Reg
        end
        JN = JN + (CNhat*Phi(:,j+1,i)-Y(i,j+1))^2; %Note that this excludes j=1, i.e. t=0, because there is no output then.
    end
    %eta(:,j,m+1)=ANhat' * eta(:,j+1,m+1) + 2*(CNhat*Phi(:,j,m+1)-test_y(j))*CNhat';
    for r=0:P
        dJN(r+M+3) = dJN(r+M+3) + eta(:,j+1,m+1)'*BNhat*SplinesP_linear(r+1,j+1);
    end
end

JN = JN + Reg(qM,c);
dJN = dJN + dReg(qM,c);

end