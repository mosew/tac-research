function eivs = evaluate_eivs(eivs,theta)
    eivs = cell2mat(cellfun(@(c) feval(c,theta),eivs,'UniformOutput',false));
end