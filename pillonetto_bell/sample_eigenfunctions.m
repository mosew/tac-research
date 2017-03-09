function sampled_eifs = sample_eigenfunctions(eifs,t,P)
    % INPUTS:
    % eifs is a 1xP cell array of function handles containing the eigenfunctions of the RK
    % t is a vector of values to sample each eigenfunction with
    % 
    % OUTPUT:
    % sampled_eifs is n x P, each column is an ef
    P=length(eifs);
    sampled_eifs = zeros(length(t),P);
    for k = 1:P
        sampled_eifs(:,k) = feval(eifs{k},t);
    end
end