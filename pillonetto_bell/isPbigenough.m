function bool = isPbigenough(fksCoarse,fksFine,T)
    % fkCoarse (resp. fkFine)
    % is a K x n array where the k^th ROW is the sum of the first P (resp. Q=floor(gamma*P)+1)
    % eigenfunctions with coeffs a^k, where k is the #iteration of the MCMC
    
    
    e = 0.02;
%    gamma = 1.25;
    
    lumC=confidence_limits(fksCoarse);
    lumF=confidence_limits(fksFine);
    
    bool = condition(e,T,lumF,lumC,1) && condition(e,T,lumF,lumC,2) && condition(e,T,lumF,lumC,3);
end