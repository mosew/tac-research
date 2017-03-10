function output = L(theta,a,efs,T)
    % I THINK IT'S BROKEN?
    % INPUTS:
    % theta is a vector of parameters
    % a,efs corresponds to the input signal to L:
    %       u(t) = a.*efs(t)
    % T is the max time.


    % efs is a cell array of P eigenfunctions
    % a is a 1-row array of P coefficients
    % th is a number.
        
    f = fk_from_a_efs(a,efs);
    
    g = @(s) exp(-theta(2)*(T-s)); %note this one is in a form ready for convolution
    
    output = integral(f*g,0,T);
end