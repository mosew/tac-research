function J = Reg_Min_Fun_L_FD(Reg_Parms)

global Training_Data
global n np1 n2 n2p2 n_state
global Regularization

global M M_INV C
global A B A_INV
global Reg_Wt_TAC Reg_Wt_BrAC Normalization

rho_M = Reg_Parms(1)^2;
rho_K = Reg_Parms(2)^2;

SS_Residual_TAC = 0;
SS_Residual_BrAC = 0;
    
t_count = Training_Data.t_count;
tau = Training_Data.tau;
u = Training_Data.u;
y = Training_Data.y;
m_discretize = Training_Data.m_discretize;
t_tau = Training_Data.t_tau;
    
A_Spline_Sim = Regularization.A_Spline_Sim;
Splines = Regularization.Splines;
A_hat = Regularization.A_hat;
B_hat = Regularization.B_hat;
Q_M = Regularization.Q_M;
Q_K = Regularization.Q_K;
    
A_Tilda = [A_Spline_Sim;rho_M*Q_M + rho_K*Q_K];
b_Tilda = [y';zeros(m_discretize+1,1)];

U_0 = zeros(m_discretize+1,1);

warning off;

OPTIONS = optimset('Display','off');
[BrAC,RESNORM,RESIDUAL,EXITFLAG] = lsqnonneg(A_Tilda,b_Tilda,U_0,OPTIONS);

warning on;

BrAC_Splines = Splines*BrAC;

x_hat = [zeros(n_state,1)];

for j = 1:t_count
    
    x_hat = [x_hat,A_hat*x_hat(:,j) + B_hat*BrAC_Splines(j,1)];
    
end

y_hat = C*x_hat;
TAS_BrAC = y_hat';

SS_Residual_TAC = SS_Residual_TAC + (norm(y-TAS_BrAC')^2);
SS_Residual_BrAC = SS_Residual_BrAC + (norm(u - BrAC_Splines')^2);
    
J = Reg_Wt_TAC*SS_Residual_TAC/Normalization + Reg_Wt_BrAC*SS_Residual_BrAC/Normalization; 

end







      
      