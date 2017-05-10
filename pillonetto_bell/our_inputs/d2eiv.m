function d2eiv = d2eiv(theta,P,T)
    % THIS IS FOR THE GREEN'S FIRST KERNEL
    
    % OUTPUT:
    % 3D matrix
    % P x nTheta x nTheta
    % Hessian of Lambda_j evaluated at theta, where the (j,s,r) entry is
    %       d^2 Lambda_j / dth_s dth_r
    
    nTheta = length(theta);
    d2eiv = zeros(P,nTheta,nTheta);
    
end