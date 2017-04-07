function infoTheta = infoTheta()
    % Gives Fisher information matrix for p(theta), i.e. the prior on theta
    % OUTPUT:
    % nTheta x nTheta
    
    % For EXAMPLE 7.2 we assume that theta(1) and theta(2) are independent,
    % so the Fisher information matrix is diagonal.
    
    % theta(1)~exponential(1)? Jeffreys prior?
    % theta(2)~N(10,2)
    
    infoTheta = [0.01,0;0,0.25];
    
end