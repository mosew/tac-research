function Vy_given_theta = Vy_given_theta(theta,T,P,rkhs_eigenfile,data_path)
    % INPUTS:
    % theta is a vector of parameters
    % T is the total length of an episode (in minutes?)
    % P is number of eigenfunctions to use
    %
    % OUTPUTS:
    % n x n matrix, where n = length(y) = T/tau
    
        
    % RKHS dependence
    [eivs,eifs] = get_kernel_eigenstuff(theta,P,T,rkhs_eigenfile);
    sampled_eifs = sample_eigenfunctions(eifs,t,P);
    
    % Data dependence
    meas_noise = Vy_given_theta_u(data_path);
    
    
    convolved_sampled_eifs = convolve_sample_eifs(sampled_eifs);
    
    
    Vy_given_theta = zeros(n+1,n+1);
    for i=1:(n+1)
        for k=1:(n+1)
            Vy_given_theta(i,k) = sum(eivs'.*convolved_sampled_eifs(:,i).*convolved_sampled_eifs(:,k)) + meas_noise(i,k);
        end
    end
    
    
end