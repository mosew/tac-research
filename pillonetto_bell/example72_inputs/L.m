function output = L(th2,a,efs,T)
    % efs is a cell array of P eigenfunctions
    % a is a 1-row array of P coefficients
    % th is a number.
    
    f = fk_from_a_efs(a,efs);
    
    g = @(s) exp(-th2*(T-s)); %note this one is in a form ready for convolution
    
    output = integral(f*g,0,T);
end