function convolved_eifs = convolve_eifs(eifs,theta,P,n,tau)
    % INPUTS:
    % eifs is a cell array of P eigenfunctions
    % theta is a parameter vector
    % P is # of eigenfunctions
    % tau is timestep = 5
    %
    % OUTPUTS:
    % convolved_eifs is a P x n array of our eigenfunctions passed to the functional L, i.e. 
    %           convolved_sampled_eifs(i,j) = L_i(theta,phi_j)

    
    % Note: right now this doesn't do anything with the sampled
    % eigenfunctions.
        
    convolved_eifs = zeros(P,n);
    
    for i = 1:n        
        for j = 1:P
            convolved_eifs(j,i) = L_i(theta,eifs{j},i,tau);
        end
    end
    
end