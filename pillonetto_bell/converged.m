function bool = converged(epsilon)

    bool = condition(epsilon,T,lumF,lumC


    bool= (k>=K);
%     if nEndingSteps<=0
%         fprintf('nEndingSteps should be > 1')
%     end
%     
%     if size(a,2)==k-1+burnoff
%         bool = 1;
%         return
%     end
%     
%     
%     if k < burnoff+nEndingSteps
%         bool=0;
%     else
%         % This condition returns true only if the amplitude guesses changed
%         % by less than epsilon the previous nEndingSteps iterations
%     bool= all(all(abs(diff(a(:,(k-burnoff-nEndingSteps+1):k-burnoff),2))<epsilon));
%     end
%     
end