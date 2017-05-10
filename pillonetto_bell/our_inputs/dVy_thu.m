function dVy_thu = dVy_thu(theta,P,n)
    % This is the derivative of the known covariance matrix
    % V(y|theta,u)
    % with respect to theta
    
    % FOR EXAMPLE 7.2 FROM PILLONETTO-BELL
    
    % The paper is not totally clear to me, but I think the v_i are iid Gaussian
    % with mean 0 and variance .05u(t) at time t
    
    % In other words, V(y|theta,u) does not depend on theta.
    
    % OUTPUT:
    % n x n x nTheta

    nTheta = length(theta);
    dVy_thu = zeros(n,n,nTheta);     
    
end