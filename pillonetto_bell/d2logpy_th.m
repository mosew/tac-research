function d2logpy_th = d2logpy_th(y_total,theta,tau,P,T,n,eivs,rkhs_eigenfile,data_path)
    
    % INPUT:
    % y_total is m x n
    % OUTPUT:
    % nTheta x nTheta, where d2logpy_th(s,r) is 
    % d (log( p(y|th) ) / d(th_s) d(th_r)
    
    
    nTheta = length(theta);
    
    Lmatrix_ = Lmatrix(theta,P,T,rkhs_eigenfile,n,tau);
    dLmatrix_ = dLmatrix(theta,P,n,Lmatrix_);

    
    m = size(y_total,1);

    s2 = d2Vy_th(theta,P,T,n,eivs,Lmatrix_,dLmatrix_);
    s1 = dVy_th(theta,P,T,n,eivs,Lmatrix_);
    s0 = Vy_th(theta,P,T,n,tau,rkhs_eigenfile,data_path);
    
    d2logpy_th = zeros(nTheta,nTheta);
    
    for  s = 1:nTheta
        for r = 1:nTheta
            d2logpy_th(s,r) = -0.5*(m*trace(s0\(s2(:,:,s,r)-(s1(:,:,s)/s0)*s1(:,:,r))) + y_total(1:m,:)*s2(:,:,s,r)*y_total(1:m,:)');
        end
    end
    
end