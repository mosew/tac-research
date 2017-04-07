function [r1_r2_h] = BrAC_Estimator_Filter_Design(BrAC,TAC)

%  System Requirements: Standard MATLAB Version 7 or Later including 
%  Optimization Toolbox

%  Input:  BrAC: First Column: Time in Hours
%                Second Column: BrAC in % Alcohol
%          TAC:  First Column: Time in Hours
%                Second Column: TAC in % Alcohol

%  Output:  r1_r2_h: Column Vector; First two entries are regularization
%                    weights, the remainder is forward convolution filter 
%                    in units of % alcohol vs time in hours for 20 hours
%                    sampled at 1 minute intervals (1/60 hours)

%  This version contains both the PC and MAC version of the code.  It is
%  controlled by a switch in the internal parameters file.  Note that if
%  being used on a MACINTOSH machine, the EXCEL file containing the TAC and
%  BRAC data must be saved in 95 mode before being used by this code

%  No additional MATLAB Toolboxes required

%  Other routines called by this program:

%  Build_B_1.m 
%  Build_C_1.m 
%  Build_K_1.m 
%  Build_L_1.m 
%  Build_M_1.m 
%  Build_Q_M.m 
%  Build_Q_K.m 
%  Reg_Min_Fun_L_FD.m
%  plot_L.m
%  plot_R.m
%  TS_J_eval_3_EXP_PARABOLIC_F_FD.m
%  TS_J_eval_3_EXP_HYPERBOLIC_F_FD.m
%  TS_J_eval_3_EXP_PARABOLIC_G_FD.m
%  TS_J_eval_3_EXP_HYPERBOLIC_G_FD.m

global Training_Data
global n np1 n2 n2p2 n_state
global Regularization

global M_1 K_1 L_1 B_1 C_1
global M M_INV C

global A B A_INV
global Reg_Wt_TAC Reg_Wt_BrAC Normalization

global CD_time_BrAC CD_data_BrAC

global DA DB



% Begin Settings

sampling_interval = 1/60;
T_IR = 20;

parms_init = zeros(1,3);

n_discretize = 32;            % n_discretize; integer recommended to be at least 32 {32}
tau_discretize = 300;         % tau_discretize {300}
parms_init(1) = 1;            % parms_init(1); first component of initial guess for parameters {1}
parms_init(2) = 1;            % parms_init(2) ); second component of initial guess for parameters {1}
parms_init(3) = 1;            % parms_init(3); third component of initial guess for parameters {1}
m_discretize_p_h = 6;         % m_discretize_p_h; recommended six (6) per hour {10}
reg_parm_M = .1;              % default regularization parameter for BrAC (reg_parm_M) {.1}
reg_parm_K =.1;               % default regularization parameter for BrAC? (reg_parm_K) {.1}
tol_fun_n = 6;                % -log(default functional convergence tolerance) {6}
tol_x_n = 6;                  % -log(default paramter convergence tolerance) {6}
Reg_Wt_TAC =.5;               % regularization weight for TAC ; should be between 0 and 1 {.5}
Reg_Wt_BrAC =.5;              % regularization weight for BrAC ; should be 1 - preceding entry {.5)
Parabolic = 1;                % model switch Parabolic 1 Hyperbolic 0
Opt_Switch = 1;               % optimization switch FMINSEARCH 0 FMINCON (Optimization Toolbox Required) 1

% End Settings

Training_Data = struct('t_count',0,'tau',0,'t_tau',[],'m_discretize',0,'y',[],'u',[]);
Regularization = struct('A_Spline_Sim',[],'Splines',[],'A_hat',[],'B_hat',[],'Q_M',[],'Q_K',[]);

BrAC_TAC_Data_BrAC_time = BrAC(:,1);
BrAC_TAC_Data_BrAC_data = BrAC(:,2);
BrAC_TAC_Data_TAC_time = TAC(:,1);
BrAC_TAC_Data_TAC_data = TAC(:,2);

CD_beg = 0;

CD_end = (max(max(BrAC_TAC_Data_BrAC_time),max(BrAC_TAC_Data_TAC_time)));

BrAC_TAC_time = CD_beg:sampling_interval:CD_end;

I_BrAC = find(((BrAC_TAC_time - (max(BrAC_TAC_Data_BrAC_time)))>= 0),1);

r_BrAC_data = max(spline(BrAC_TAC_Data_BrAC_time,BrAC_TAC_Data_BrAC_data,...
    BrAC_TAC_time(1:I_BrAC)),0);
n_BrAC = length(BrAC_TAC_time)-I_BrAC;

if (n_BrAC ~=0)
    
    r_BrAC_data = [r_BrAC_data,r_BrAC_data(end) - (1:n_BrAC)*r_BrAC_data(end)/n_BrAC]; 
    
end
    
I_TAC = find(((BrAC_TAC_time - (max(BrAC_TAC_Data_TAC_time)))>= 0),1);

r_TAC_data = max(1000*spline(BrAC_TAC_Data_TAC_time,BrAC_TAC_Data_TAC_data,...
   BrAC_TAC_time(1:I_TAC)),0);
n_TAC = length(BrAC_TAC_time)-I_TAC;

if (n_TAC ~= 0)
    
    r_TAC_data = [r_TAC_data,r_TAC_data(end) - (1:n_TAC)*r_TAC_data(end)/n_TAC]; 
    
end

CD_time_BrAC = BrAC_TAC_time;
CD_time_TAC = BrAC_TAC_time;

CD_data_TAC = r_TAC_data;

CD_data_BrAC = 1000*r_BrAC_data;

t_BrAC = CD_time_BrAC';
data_BrAC = CD_data_BrAC;

t_TAC = CD_time_TAC';
data_TAC = CD_data_TAC;

if (t_TAC(end) > t_BrAC(end))
    
    tau_bar_BrAC = mean(diff(t_BrAC));
    n_BrAC = ceil((t_TAC(end) - t_BrAC(end))/tau_bar_BrAC);
    t_BrAC = [t_BrAC,t_BrAC(end) + (1:n_BrAC)*((t_TAC(end) - t_BrAC(end))/n_BrAC)];
    
    data_BrAC = [data_BrAC,data_BrAC(end) - (1:n_BrAC)*data_BrAC(end)/n_BrAC];
     
end

if (t_TAC(end) < t_BrAC(end))
    
    tau_bar_TAC = mean(diff(t_TAC));
    n_TAC = ceil((t_BrAC(end) - t_TAC(end))/tau_bar_TAC);
    t_TAC = [t_TAC,t_TAC(end) + (1:n_TAC)*((t_BrAC(end) - t_TAC(end))/n_TAC)];
    
    data_TAC = [data_TAC,data_TAC(end) - (1:n_TAC)*data_TAC(end)/n_TAC];
    
end

if (tau_discretize ~= 0)
    
    Training_Data.t_count = tau_discretize;
    
else
    
    Training_Data.t_count = length(t_TAC) - 1;
    
end

Training_Data.tau = (t_TAC(end) - t_TAC(1))/Training_Data.t_count;
Training_Data.t_tau = (0:Training_Data.t_count)*Training_Data.tau;

y_spline_coef = spline(t_TAC,data_TAC');
Training_Data.y = ppval(y_spline_coef,Training_Data.t_tau);

Training_Data.u = interp1(t_BrAC,data_BrAC',Training_Data.t_tau,'linear','extrap');

Training_Data.m_discretize = ceil(m_discretize_p_h * Training_Data.t_tau(end));

tau = Training_Data.tau;

t_tau_IR = (0:ceil(T_IR/sampling_interval))*sampling_interval;

n = n_discretize;
n2 = 2*n;
np1 = n + 1;
n2p2 = n2 + 2;

if (Parabolic)
    
    n_state = np1;
    
else
    
    n_state = n2p2;
    
end

K_1 = Build_K_1(n);
L_1 = Build_L_1(n);
B_1 = Build_B_1(n);

if (Parabolic)
    
    M = Build_M_1(n); 
    M_INV = inv(M);
    
    C = Build_C_1(n);
    
else
    
    M_1 = Build_M_1(n); 
    C_1 = Build_C_1(n);
    
    M = [M_1,zeros(np1,np1);zeros(np1,np1),M_1];
    M_1_INV = inv(M_1);
    M_INV = [M_1_INV,zeros(np1,np1);zeros(np1,np1),M_1_INV];

    C = [C_1,zeros(1,np1)];
    
end

D = [0];
          
if (Parabolic)
           
   parms_init = parms_init(1:2);
             
   if (Opt_Switch)
               
      DA = -M_INV*K_1;
      DB = M_INV*B_1;
  
      Aeq = [];
      Beq = [];
      LB = zeros(size(parms_init));
      UB = inf*ones(size(parms_init));
      NONLCON = [];
      Aineq = [];
      Bineq = [];

      OPTIONS = optimset('Display','off','GradObj','on','MaxIter',1000,'MaxFunEvals',3000,'TolFun',10^-tol_fun_n,...
                    'TolX',10^-tol_x_n);
      warning off
      [parms_star,FVAL,EXITFLAG] = fmincon('TS_J_eval_3_EXP_PARABOLIC_G_FD',parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,...
                    NONLCON,OPTIONS);
      warning on
      alpha = parms_star(1);
      beta  = parms_star(2);
                
   else
       
      OPTIONS = optimset('Display','off','MaxIter',1000,'MaxFunEvals',3000,'TolFun',10^-tol_fun_n,'TolX',10^-tol_x_n);
      warning off
      [parms_star,FVAL,EXITFLAG] = fminsearch('TS_J_eval_3_EXP_PARABOLIC_F_FD',parms_init,OPTIONS);
                warning on
        
      alpha = parms_star(1)^2;
      beta  = parms_star(2)^2;
               
   end

   K = -alpha*K_1 - L_1;
   F = beta*B_1;
                        
else
           
   if (Opt_Switch)
               
      DA = zeros(n2p2,n2p2,3);
      DB = zeros(n2p2,3);
      
      DA(np1+1:n2p2,1:np1,1) = -M_1_INV*K_1;
      DA(np1+1:n2p2,np1+1:n2p2,2) = -eye(np1,np1);
                
      DB(np1+1:n2p2,3) = M_1_INV*B_1;
   
      Aeq = [];
      Beq = [];
      LB = zeros(size(parms_init));
      UB = inf*ones(size(parms_init));
      NONLCON = [];
      Aineq = [];
      Bineq = [];
                
               
      OPTIONS = optimset('Display','off','GradObj','on','MaxIter',1000,'MaxFunEvals',3000,...
                    'TolFun',10^-tol_fun_n,'TolX',10^-tol_x_n);
      warning off
      [parms_star,FVAL,EXITFLAG] = fmincon('TS_J_eval_3_EXP_HYPERBOLIC_G_FD',parms_init,Aineq,...
                    Bineq,Aeq,Beq,LB,UB,NONLCON,OPTIONS);
      warning on
                
      alpha = parms_star(1);
      beta  = parms_star(2);
      gamma = parms_star(3);
                
   else
               
      OPTIONS = optimset('Display','off','MaxIter',1000,'MaxFunEvals',3000,'TolFun',10^-tol_fun_n,'TolX',10^-tol_x_n);
      warning off
      [parms_star,FVAL,EXITFLAG] = fminsearch('TS_J_eval_3_EXP_HYPERBOLIC_F_FD',parms_init,OPTIONS);
                warning on
        
      alpha = parms_star(1)^2;
      beta  = parms_star(2)^2;
      gamma = parms_star(3)^2;
               
   end
                
   K = [zeros(np1,np1),M_1;-alpha*K_1-L_1,-beta*M_1];
   F = [zeros(np1,1); gamma*B_1];
                    
end
       
A = M_INV*K;
A_INV = inv(A);

B = M_INV*F;
       
A_hat = expm(tau*A);
B_hat = A_INV*(A_hat - eye(n_state))*B;
               
SYS = ss(A,B,C,[0]);
            
[Impulse_Response_y,Impulse_Response_t] = impulse(SYS,t_tau_IR);

t_count = Training_Data.t_count;
tau = Training_Data.tau;
u = Training_Data.u;
y = Training_Data.y;
m_discretize = Training_Data.m_discretize;
t_tau = Training_Data.t_tau;

A_hat = expm(tau*A);
B_hat = A_INV*(A_hat - eye(n_state))*B;

Spline_nodes = (0:m_discretize)*t_tau(end)/m_discretize;
Spline_Values = eye(m_discretize+1);
Splines = interp1(Spline_nodes,Spline_Values,t_tau);

A_Spline_Sim = [];

for i=1:m_discretize+1,
    
    x_hat = [zeros(n_state,1)];
    
    for j = 1:t_count
        
        x_hat = [x_hat,A_hat*x_hat(:,j) + B_hat*Splines(j,i)];
        
    end
    
    y_hat = C*x_hat;
    
    A_Spline_Sim = [A_Spline_Sim,y_hat'];
    
end

Q_M = t_tau(end)*Build_Q_M(m_discretize);
Q_K = (1/t_tau(end))*Build_Q_K(m_discretize);

Regularization.A_Spline_Sim = A_Spline_Sim;
Regularization.Splines = Splines;
Regularization.A_hat = A_hat;
Regularization.B_hat = B_hat;
Regularization.Q_M = Q_M;
Regularization.Q_K = Q_K;

reg_parm_M_Inv = reg_parm_M;
reg_parm_K_Inv = reg_parm_K;

A_Tilda = [A_Spline_Sim;reg_parm_M*Q_M + reg_parm_K*Q_K];
b_Tilda = [y';zeros(m_discretize+1,1)];

U_0 = zeros(m_discretize+1,1);

OPTIONS = optimset('Display','off');
warning off
[BrAC,RESNORM,RESIDUAL,EXITFLAG] = lsqnonneg(A_Tilda,b_Tilda,U_0,OPTIONS);
warning on

BrAC_Splines = Splines*BrAC;

x_hat = [zeros(n_state,1)];

for j = 1:t_count
    
    x_hat = [x_hat,A_hat*x_hat(:,j) + B_hat*BrAC_Splines(j,1)];
    
end

y_hat = C*x_hat;

TAS_BrAC = y_hat';

SS_Residual_TAC = (norm(y-TAS_BrAC')^2);
SS_Residual_BrAC_Int = (norm(u - BrAC_Splines')^2);

Normalization = SS_Residual_TAC + SS_Residual_BrAC_Int;
         
reg_parms_init = [sqrt(reg_parm_M),sqrt(reg_parm_K)];
    
OPTIONS = optimset('Display','off','MaxIter',1000,'MaxFunEvals',3000,'TolFun',10^-8,'TolX',10^-6);
warning off
[reg_parms_star,FVAL,EXITFLAG] = fminsearch('Reg_Min_Fun_L_FD',reg_parms_init,OPTIONS);
warning on
      
reg_parm_M_Inv = reg_parms_star(1)^2;
reg_parm_K_Inv = reg_parms_star(2)^2;

r1_r2_h = [reg_parm_M_Inv,reg_parm_K_Inv,Impulse_Response_y'];
    
end
 