function dVy_thu = dVy_thu(theta,P,n)
    % This is the derivative of the known covariance matrix
    % V(y|theta,u)
    % with respect to theta
    
    % FOR ROSEN TRANSDERMAL ALCOHOL MODEL
        
    % In other words, V(y|theta,u) does not depend on theta.
    
    % OUTPUT:
    % n x n x nTheta

    nTheta = length(theta);
    dVy_thu = zeros(n,n,nTheta);     
    
end