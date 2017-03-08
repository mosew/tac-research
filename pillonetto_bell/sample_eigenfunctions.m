function sampled_eifs = sample_eigenfunctions(eifs,t,P)
    % Output is M x P, each column is an ef
    sampled_eifs = zeros(length(t),P);
    for k = 1:P
        sampled_eifs(:,k) = feval(eifs{k},t);
    end
end