% main routine

clear
JN_and_dJN_globals


% Test parameters
qtest = ones(1,M+2);
utest = zeros(1,n);

parms_init = [qtest,utest];



% Options for constrained minimization
Aineq = [];
Bineq = [];
LB = [eps*ones(size(qtest)),zeros(size(utest))];
UB = inf*ones(size(parms_init));
NONLCON = [];
Aeq = [];
Beq = [];
OPTIONS = optimset('Display','off','GradObj','on','MaxIter',5000,'MaxFunEvals',3000,'TolFun',10^-9,...
            'TolX',10^-9);


% Constrained minimization        
[parms_star,FVAL,EXITFLAG,OUTPUT] = fmincon('JN_and_dJN',parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,...
            NONLCON,OPTIONS);

        
        
% Collect optimal parameters and deconvolved signal
q2_star = parms_star(1)
q1M_star = parms_star(2:M+2)
u_star  = parms_star(M+3:end);





% Plot performance
figure
hold on
%plot(u_star)


% % Post-processing: Moving Average
% b=14;
% u_star_ma = u_star;
% u_s = [zeros(1,b/2),u_star,zeros(1,b/2)];
% for i = b/2+1 : length(u_s)-b
%     u_star_ma(i-b/2) = mean(u_s(i-b/2:i+b/2));
% end
% plot(u_star_ma)


% % Post-processing: Savitsky-Golay
% u_star_sg = max(sgolayfilt(u_star',2,13),0);
% plot(u_star_sg)

xlim([0,length(test_y)])
plot(test_u,'.')
%plot(test_y)


% Feed optimal parameters and input forward through system to check it
figure
forward_system
