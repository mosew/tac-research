function dLmatrix2 = dLmatrix2(theta,P,n,Lmatrix)
    % OUTPUT:
    % P x n x nTheta, where dLmatrix(j,i,r) = dL_i(th,phi_j) / dth_r
    
    % ROSEN TRANSDERMAL ALCOHOL MODEL
    
    nTheta = length(theta);    
    
    load('operators.mat','dAN_dq1');
    
    dLmatrix2 = zeros(P,n,nTheta);
%     dLmatrix(:,:,1) = dAN_dq1*Lmatrix;
    dLmatrix2(:,:,2) = Lmatrix/theta(2);
    
    
    
end