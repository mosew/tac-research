function p_theta = p_theta(i,theta)
    % this is for example 7.2 from the Pillonetto-Bell paper
    
    switch i
        case 1
            % p_theta = p_y_given_theta(y,theta,tau,T,P,rkhs_eigenfile,data_path);
            % flat prior

            if theta(1)<=0 || theta(1)>=15
                p_theta = 0;
            else
                p_theta=1/15;
            end
        case 2
            if theta(2)<=0
                p_theta = 0;
            else
                p_theta=normpdf(theta(2),10,2);
            end
    end
end