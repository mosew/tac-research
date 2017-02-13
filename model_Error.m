function ss = model_Error(q2,qM,c,test_u,V,P,T)
    global CNhat
    [eiv,eif]=rkhs_eigen(theta,P,T);
    b = a./sqrt(eiv);
    
    C = zeros(n+1,P);
    
    for i = 1:length(eif)
        Phi = forward_system(q2,qM,eif{i}(0:1/n:1));
        for j = 0:n
            C(j+1,i) = CNhat*Phi(:,j+1,m+1);
        end
    end
    
    e = C*b-test_u;
    ss = 0.5*sum(b.^2) + 0.5*(e'/V)*e;
    
end