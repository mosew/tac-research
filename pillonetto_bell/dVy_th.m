function dVy_th = dVy_th(theta,P,T,n,eivs,Lmatrix)
    % OUTPUT:
    % 3D matrix
    % n x n x nTheta
    % 
    % dVy_th(i,k,r) =  ( dV(y|th)/dth_r )_ik
    
    eivs = evaluate_eivs(eivs,theta);
    
    nTheta = length(theta);
    
    deiv_ = deiv(theta,P,T);
        
    dLmatrix_=dLmatrix(theta,P,n,Lmatrix);
    dVy_thu_ = dVy_thu(theta,P,n);
    
    dVy_th = zeros(n,n,nTheta);
    
    for i = 1:n
        for k = 1:n
            for r = 1:nTheta
                dVy_th(i,k,r) = sum( deiv_(:,r).*Lmatrix(:,i).*Lmatrix(:,k) + eivs.*(dLmatrix_(:,i,r).*Lmatrix(:,k)+Lmatrix(:,i).*dLmatrix_(:,k,r))) ;
            end
        end
    end
    
    dVy_th = dVy_th + dVy_thu(theta,P,n);

end