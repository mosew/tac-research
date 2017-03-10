function convolved_sampled_eifs = convolve_sampled_eifs(sampled_eifs,P)

    convolved_sampled_eifs = zeros(size(sampled_eifs));
    for i = 1:P
        convolved_sampled_eifs(i,:) = L(theta,a,efs,T);
    end
    
end