function lo_up_mid = confidence_limits(fks)
    % fks is a 1 x K array of function handles
    % where the k^th column represents the function handle corresponding to
    % the k^th run (the k^th amplitude draws) from the MCMC
        
    d=@(t) cell2mat(cellfun(@(c) feval(c,t),fks','UniformOutput',false));    
    
    fL = @(t) quantile(feval(d,t),0.025);
    fU = @(t) quantile(feval(d,t),0.975);
    fM = @(t) mean(feval(d,t));
    lo_up_mid={fL,fU,fM};
end