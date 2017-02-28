function sampled_eifs = sample_eigenfunctions(eifs,t,P)
    % Output is n x P, each row is an ef
    sampled_eifs = zeros(P,length(t));
    for k = 1:P
        sampled_eifs(k,:) = feval(eifs{k},t);
    end
end