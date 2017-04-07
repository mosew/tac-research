function bool = converged(a,k,nEndingSteps,burnout,epsilon)
    if nEndingSteps<=0
        fprintf('nEndingSteps should be > 1')
    end
    
    if size(a,2)==k-1+burnout
        bool = 1;
        return
    end
    
    
    if k < burnout+nEndingSteps
        bool=0;
    else
        % This condition returns true only if the amplitude guesses changed
        % by less than epsilon the previous nEndingSteps iterations
    bool= all(all(abs(diff(a(:,(k-burnout-nEndingSteps+1):k-burnout),2))<epsilon));
    end
    
end