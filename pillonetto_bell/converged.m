function bool = converged(a,k,nEndingSteps,epsilon)
    if nEndingSteps<=0
        fprintf('nEndingSteps should be > 1')
    end
    
    if k < nEndingSteps
        bool=0;
    else
        % This condition returns true only if the amplitude guesses changed
        % by less than epsilon the previous nEndingSteps iterations
        
    bool= all(all( abs(diff(a(:,(k-nEndingSteps+1):k),2))<epsilon));
    end
end