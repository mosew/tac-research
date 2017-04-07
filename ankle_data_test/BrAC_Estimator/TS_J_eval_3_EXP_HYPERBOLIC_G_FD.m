function [J,grad_J] = TS_J_eval_3_EXP_HYPERBOLIC_G_FD(parms)

global Training_Data
global n np1 n2 n2p2 n_state

global M_1 K_1 L_1 B_1 C_1
global M M_INV C

global DA DB

alpha = parms(1);
beta  = parms(2);
gamma = parms(3);

K = [zeros(np1,np1),M_1;-alpha*K_1 - L_1,-beta*M_1];
F = [zeros(np1,1); gamma*B_1];

A = M_INV*K;
A_INV = inv(A);

B = M_INV*F;
    
t_count = Training_Data.t_count;
tau = Training_Data.tau;
u = Training_Data.u;
y = Training_Data.y;

A_EXP_alpha = expm([tau*A,tau*DA(:,:,1);zeros(n2p2,n2p2),tau*A]);
A_EXP_beta = expm([tau*A,tau*DA(:,:,2);zeros(n2p2,n2p2),tau*A]);

A_hat = A_EXP_alpha(1:n2p2,1:n2p2);
B_hat = A_INV*(A_hat - eye(n_state))*B;

DA_hat_alpha = A_EXP_alpha(1:n2p2,n2p2+1:2*n2p2);
DA_hat_beta = A_EXP_beta(1:n2p2,n2p2+1:2*n2p2);

DB_hat_alpha = -A_INV*(DA(:,:,1)*A_INV*(A_hat - eye(n_state)) - DA_hat_alpha)*B;
DB_hat_beta = -A_INV*(DA(:,:,2)*A_INV*(A_hat - eye(n_state)) - DA_hat_beta)*B;
DB_hat_gamma = A_INV*(A_hat - eye(n_state))*DB(:,3);

x = [zeros(n_state,1)];

for j = 1:t_count
    
    x = [x,A_hat*x(:,j) + B_hat*u(:,j)];
    
end

e = C*x - y;

J = trace(e'*e);

grad_J = zeros(1,3);

z = [zeros(n_state,1)];

for j = t_count:-1:1
    
    z = [z,A_hat'*z(:,t_count-j+1) + 2*C'*e(:,j+1)];
    
    grad_J(1,1) = grad_J(1,1) + z(:,t_count-j+2)'*(DA_hat_alpha*x(:,j) + DB_hat_alpha*u(:,j));
    grad_J(1,2) = grad_J(1,2) + z(:,t_count-j+2)'*(DA_hat_beta*x(:,j) + DB_hat_beta*u(:,j));
    grad_J(1,3) = grad_J(1,3) + z(:,t_count-j+2)'*DB_hat_gamma*u(:,j);
    
end
     
end



