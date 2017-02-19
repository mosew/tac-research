% main routine

JN_and_dJN_globals
global SplinesP_linear
global M P


% Do we want a constant diffusivity parameter?
if M==0
    M=1;
    JN_and_dJN_globals
    Aeq = [0,1,-1,zeros(1,P+1)];
    Beq = 0;
else
    Aeq=[];
    Beq=[];
end

% Test parameters
qtest = ones(1,M+2);
ctest = zeros(1,P+1);

parms_init = [qtest,ctest];



% Options for constrained minimization
Aineq = [];
Bineq = [];
LB = [eps*ones(size(qtest)),zeros(size(ctest))];
UB = [inf*ones(size(qtest)),inf*ones(size(ctest))];
NONLCON = [];
OPTIONS = optimset('Display','off','GradObj','on','MaxIter',5000,'MaxFunEvals',3000,'TolFun',10^-7,...
            'TolX',10^-7);

        
        
badscale=0;
lastwarn('')
% Constrained minimization        
[parms_star,FVAL,EXITFLAG,OUTPUT] = fmincon('JN_and_dJN_uspline',parms_init,Aineq,Bineq,Aeq,Beq,LB,UB,...
            NONLCON,OPTIONS);
if lastwarn
    badscale=1;
    lastwarn('')
end
        
        
% Collect optimal parameters and deconvolved signal
q2_star = parms_star(1);
q1M_star = parms_star(2:M+2);
u_star  = parms_star(M+3:end)*SplinesP_linear;




% % Plot performance
% figure
% hold on
% plot(u_star)
% plot(u_total(test,:))
% 
% 
% % % Post-processing: Moving Average
% % b=8;
% % u_star_ma = u_star;
% % u_s = [zeros(1,b/2),u_star,zeros(1,b/2)];
% % for i = b/2+1 : length(u_s)-b
% %     u_star_ma(i-b/2) = mean(u_s(i-b/2:i+b/2));
% % end
% % plot(u_star_ma)
% 
% 
% % % Post-processing: Savitsky-Golay
% % u_star_sg = max(sgolayfilt(u_star',2,13),0);
% % plot(u_star_sg)
% 
% xlim([0,length(test_u)])
% plot(test_u)
% %plot(test_y)


% % Feed optimal parameters and input forward through system to check it
% figure
% Phi = forward_system(q2_star,q1M_star,u_star);
% y_out = zeros(1,n);
% for j=1:n
%     y_out(j) = CNhat * Phi(:,j,1);
% end
% plot(y_out)
% hold on
% plot(test_y)