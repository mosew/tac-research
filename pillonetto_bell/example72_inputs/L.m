function output = L(th2,a,efs,T)
    % efs is a cell array of P eigenfunctions
    % a is a 1-row array of P coefficients
    % th is a number.
    
    f = @(s) sum(a.*cellfun(@(c) feval(c,s),efs));
    g = @(s) exp(-th2*(T-s));
    
    output = integral(f*g,0,T);
end