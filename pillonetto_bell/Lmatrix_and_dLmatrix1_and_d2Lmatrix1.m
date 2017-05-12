function [Lmatrix_,dLmatrix1_,d2Lmatrix1_]=Lmatrix_and_dLmatrix1_and_d2Lmatrix1(N,theta,P,T,rkhs_eigenfile,n,tau)
    % OUTPUT:
    % P x n
    % 
    % Lmatrix_(j,i) = L_i(theta,phi_j), starting everything at t=0.
    % dLmatrix1_(j,i) = dL_i/dq1(theta,phi_j), starting everything at t=0.
    
    [~,eifs] = get_kernel_eigenstuff(P,T,rkhs_eigenfile);
    
    Lmatrix_ = zeros(P,n);
    dLmatrix1_ = zeros(P,n);
    d2Lmatrix1_ = zeros(P,n);
    load('operators.mat','dAN_dq1','ANhat','BNhat','d2AN_dq1');
    kernel = @(s) mpower(ANhat,fix(s/tau)-1) * BNhat;
    
    d2m = d2AN_dq1+dAN_dq1^2;
    
    for i = 1:n
        i
        for j = 1:P
            f=eifs{j};
            fxg = @(s) kernel(i*tau-s).*feval(f,s);
            output = integral(fxg,0,i*tau,'ArrayValued',true);
            
            Lmatrix_(j,i) = output(1);
            
            doutput = (i*tau-1)*dAN_dq1*output;
            dLmatrix1_(j,i) = doutput(1);
            d2output = d2m*output;
            d2Lmatrix1_(j,i) = d2output(1);
        end
    end
    
end