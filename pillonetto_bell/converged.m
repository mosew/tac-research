function bool = converged(a,k,nEndingSteps,epsilon)
    if nEndingSteps<=0
        fprintf('nEndingSteps should be > 1')
    end
    
    if k < nEndingSteps
        bool=0;
    else
        % This condition returns true only if the amplitude guesses changed
        % by less than epsilon the previous nEndingSteps iterations
        bool=max(abs(diff(a(:,k:-1:(k-nEndingSteps)),2)))<epsilon;
    end
end