function p_theta = p_theta(i,theta)
    % this is for example 7.2 from the Pillonetto-Bell paper
    
    switch i
        case 1
            % p_theta = p_y_given_theta(y,theta,tau,T,P,rkhs_eigenfile,data_path);
            % flat/logarithmic prior
            % In this case it's OK because for MCMC we only need to find a
            % function proportional to a pdf.
            if theta(1)<=0
                p_theta = 0;
            else
                p_theta=log(theta(1));
            end
        case 2
            p_theta=normpdf(theta(2),10,2);
    end
end