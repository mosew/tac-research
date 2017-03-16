function d2eiv = d2eiv(theta,P,T)
    % THIS IS FOR EXAMPLE 7.2 in PB
    
    % OUTPUT:
    % 3D matrix
    %
    % Hessian of Lambda_j evaluated at theta, where the (j,s,r) entry is
    %       d^2 Lambda_j / dth_s dth_r
    
    nTheta = length(theta);
    
    d2eiv = zeros(P,nTheta,nTheta);
    
end