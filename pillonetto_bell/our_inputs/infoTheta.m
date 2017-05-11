function infoTheta = infoTheta()
    % Gives Fisher information matrix for p(theta), i.e. the prior on theta
    % OUTPUT:
    % nTheta x nTheta
    
    % For ROSEN TRANSDERMAL ALCOHOL MODEL we assume that theta(1) and theta(2) are independent,
    % so the Fisher information matrix is diagonal.
    
    % theta(1)~N(4.4e-3,4e-4)
    % theta(2)~N(1.23,0.1)
    
    infoTheta = [1.6e-7,0;0,1e-2];
    
end