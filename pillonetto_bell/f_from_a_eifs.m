function f = f_from_a_eifs(a,eifs)

    % INPUTS:
    % a is a 1xP vector of coefficients
    % eifs is a Px1 cell array of eigenfunction handles
    %
    % OUTPUTS:
    % f is a function handle representing 
    %       f(t) = sum_{j=1}^P   a_j * phi_j(t)
    f = @(s) a*cell2mat(cellfun(@(c) feval(c,s),eifs,'UniformOutput',false));
end