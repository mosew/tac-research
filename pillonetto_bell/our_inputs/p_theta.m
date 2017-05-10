function p_theta = p_theta(theta)
    % this is for example 7.2 from the Pillonetto-Bell paper
    % Returns number, joint prior probability for theta evaluated at theta
    
    p_theta = (theta(1)>0)*(theta(2)>0)*normpdf(theta(2),10,2)/(1-normcdf(-5,0,1));%*exp(-1000*theta(1));

end