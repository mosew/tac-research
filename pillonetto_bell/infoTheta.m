function infoTheta = infoTheta()
    % Gives Fisher information matrix for p(theta), i.e. the prior on theta
    % OUTPUT:
    % nTheta x nTheta
    
    % For EXAMPLE 7.2 we assume that theta(1) and theta(2) are independent,
    % so the Fisher information matrix is diagonal.
    
    % theta(1)~exponential(1)?
    % theta(2)~N(10,2)
    
    infoTheta = [1,0;0,1/4];
    
end