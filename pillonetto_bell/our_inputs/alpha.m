function a = alpha(js)
    % For 2nd Green's kernel calculation
    % Gives array of alpha_j for prescribed j values in vector js
    fun=@(x) 1/cosh(x)+cos(x);
    a = zeros(1,length(js));
    for j=1:length(js)
        a(j) = fzero(fun,(js(j)-1/2)*pi);
    end
end