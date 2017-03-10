function output = L_current_a(theta,a,efs,T)
    % I THINK IT'S BROKEN?
    % INPUTS:
    % theta is a vector of parameters
    %   This is for example 7.2
    %       theta(1) and theta(2)
    % (a,efs) are 1xP (number,cell) arrays corresponds to the input signal to L:
    %       u(t) = a.*efs(t)
    % T is the max time.

    % efs is a cell array of P eigenfunctions
    % a is a 1-row array of P coefficients
        
    f = fk_from_a_efs(a,efs);
    
    g = @(s) exp(-theta(2)*(T-s)); %note this one is in a form ready for convolution
    
    output = integral(f*g,0,T);
end