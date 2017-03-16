function dLmatrix = dLmatrix(theta,P,n,Lmatrix)
    % OUTPUT:
    % P x n x nTheta, where dLmatrix(j,i,r) = dL_i(th,phi_j) / dth_r
    
    % EXAMPLE 7.2 FROM PILLONETTO-BELL
    
    nTheta = length(theta);    
    
    dLmatrix = zeros(P,n,nTheta);
    dLmatrix(:,:,2) = -2*theta(2)*Lmatrix;
    
    
    
end