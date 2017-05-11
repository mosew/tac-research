function d2Lmatrix = d2Lmatrix(theta,P,n,Lmatrix)
    % OUTPUT:
    % P x n x nTheta x nTheta matrix, where
    %
    % d2Lmatrix(j,i,s,r) = d^2 L_i(th,phi_j) / dth_s dth_r
    
    % FOR ROSEN'S TRANSDERMAL ALCOHOL MODEL
    
    load('operators.mat','dAN_dq1','d2AN_dq1');
    
    nTheta = length(theta);
    
    d2Lmatrix = zeros(P,n,nTheta,nTheta);
    
    d2Lmatrix(:,:,1,1) = (d2AN_dq1+dAN_dq1^2)*Lmatrix;
    d2Lmatrix(:,:,1,2) = dAN_dq1/theta(2)*Lmatrix;
    d2Lmatrix(:,:,2,1) = d2Lmatrix(:,:,1,2);
    


end