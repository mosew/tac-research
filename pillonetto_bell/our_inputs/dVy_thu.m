function dVy_thu = dVy_thu(theta,tau,P,n)
    % This is the derivative of the known covariance matrix
    % V(y|theta,u)
    % with respect to theta
    
    % FOR ROSEN TRANSDERMAL ALCOHOL MODEL
        
    % V(y|theta,u) depends exponentially on theta(1) and linearly on theta(2)
    
    % OUTPUT:
    % n x n x nTheta

    nTheta = length(theta);
    dVy_thu = zeros(n,n,nTheta);
    load('operators.mat','tau')
    dVy_thu(:,:,1) = 0:(n-1)*tau
    
end