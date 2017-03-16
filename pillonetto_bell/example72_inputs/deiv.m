function deiv = deiv(theta,P,T)
    % This is for Example 7.2 of Pillonetto-Bell
    
    % OUTPUT:
    % P x nTheta array
    %   each row/column (j,i) is
    %       dLambda_j(theta)/d(theta_i)
    %   where Lambda_j is the j^th eigenvalue of the reproducing kernel.
    
    nTheta = length(theta);
    deiv = zeros(P,nTheta);
    
    js = 1:P;
    
    deiv(:,1) = T^2 ./ ( (js-1)*pi+pi/2).^2;
    deiv(:,2) = 0;
end