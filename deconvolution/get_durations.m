function BrAC = trim_BrAC(BrAC)
    % BrAC should be an m x n matrix, with m episodes of length n each    
    nothing = zeros(m,1);
    while all(BrAC(:,end)==nothing)
        BrAC = BrAC(:,1:(end-1));
    end
end