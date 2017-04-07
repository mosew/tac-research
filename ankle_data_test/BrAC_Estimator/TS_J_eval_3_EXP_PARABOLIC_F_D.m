function J = TS_J_eval_3_EXP_PARABOLIC_F_D(parms)

global tau t_count t_tau
global u y
global n np1 n2 n2p2 n_state

global M_1 K_1 L_1 B_1 C_1
global M M_INV C

alpha = parms(1)^2;
beta  = parms(2)^2;

K = -alpha*K_1-L_1;
Gamma_Plus = (beta/(1+alpha))*(B_1 + alpha);

A_1 = M_INV*K;
A_hat_1 = expm(tau*A_1);

A_hat = [A_hat_1,(1/(6*n))*A_hat_1*M_INV(:,n);zeros(1,n),0];
B_hat = (eye(n_state) - A_hat)*Gamma_Plus;

x = zeros(n_state,1);

for j = 1:t_count
   
   x = [x,A_hat*x(:,j) + B_hat*u(:,j)];
   
end

e = C*x - y;

J = trace(e'*e);







