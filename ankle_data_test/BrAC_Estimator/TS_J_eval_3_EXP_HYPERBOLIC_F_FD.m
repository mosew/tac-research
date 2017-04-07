function J = TS_J_eval_3_EXP_HYPERBOLIC_F_FD(parms)

global Training_Data
global n np1 n2 n2p2 n_state 

global M_1 K_1 L_1 B_1 C_1
global M M_INV C

alpha = parms(1)^2;
beta  = parms(2)^2;
gamma = parms(3)^2;

K = [zeros(np1,np1),M_1;-alpha*K_1 - L_1,-beta*M_1];
F = [zeros(np1,1); gamma*B_1];

A = M_INV*K;
A_INV = inv(A);

B = M_INV*F;

t_count = Training_Data.t_count;
tau = Training_Data.tau;
u = Training_Data.u;
y = Training_Data.y;

A_hat = expm(tau*A);
B_hat = A_INV*(A_hat - eye(n_state))*B;

x = [zeros(n_state,1)];

for j = 1:t_count
    
    x = [x,A_hat*x(:,j) + B_hat*u(:,j)];
    
end

e = C*x - y;

J = trace(e'*e);

end









