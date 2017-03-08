function bool = isPbigenough(fkCoarse,fkFine,T)
    % fkCoarse and fkFine are function handles, where
    % fkCoarse is the sum of the first, say, P eigenfunctions with coeffs a^k, where k is the #iteration of the MCMC
    % fkFine is the sum of the first int(gamma*P)+1 eigenfunctions with coeffs a^k, where k is the #iteration of the MCMC

    e = 0.02;
%    gamma = 1.25;
    
    lumC=confidence_limits(fkCoarse);
    lumF=confidence_limits(fkFine);        
    
    bool = condition(e,T,lumF,lumC,1) && condition(e,T,lumF,lumC,2) && condition(e,T,lumF,lumC,3);
end