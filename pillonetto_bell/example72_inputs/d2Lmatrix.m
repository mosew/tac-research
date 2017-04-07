function d2Lmatrix = d2Lmatrix(theta,P,n,Lmatrix)
    % OUTPUT:
    % P x n x nTheta x nTheta matrix, where
    %
    % d2Lmatrix(j,i,s,r) = d^2 L_i(th,phi_j) / dth_s dth_r
    
    % FOR EXAMPLE 7.2 FROM PILLONETTO-BELL
    
    nTheta = length(theta);
    
    d2Lmatrix = zeros(P,n,nTheta,nTheta);
    
    d2Lmatrix(:,:,1,2) = -2*theta(2)/theta(1)*Lmatrix;
    d2Lmatrix(:,:,2,1) = d2Lmatrix(:,:,1,2);
    
    d2Lmatrix(:,:,2,2) = (4*theta(2)^2-2)*Lmatrix;


end