% Forward System

qu = [q2_star,q1M_star,u_star];
global M n
global test_y
global CNhat

% Process inputs
q2 = qu(1);
%q3 = qM_and_u(2);
% qM is a vector of linear spline coefficients.
qM = qu(2:M+2);
u = qu((M+3):end);
assert(length(u)==n);


Kq = build_Kq(qM);
AN = build_AN(Kq,q2);
BN = build_BN(q2);
[ANhat,dANhat_dqM] = build_expm_stuff(AN);
BNhat = build_BNhat(AN,ANhat,BN);
dBNhat_dqM = build_dBNhat_dqM(AN,ANhat,dANhat_dqM,BN);
Phi = build_Phi(ANhat,BNhat,u);

y_out = zeros(1,n);
for j=1:n
    y_out(j) = CNhat * Phi(:,j,1);
end
plot(y_out)
hold on
plot(test_y)