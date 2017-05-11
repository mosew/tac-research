function p_theta = p_theta(theta)
    % this is for ROSEN TRANSDERMAL ALCOHOL MODEL
    % Returns number, joint prior probability for theta evaluated at theta
    
    p_theta = theta(1)>0 * theta(2)>0 * normpdf(theta(1),4.4e-3,4e-4)*normpdf(theta(2),1.23,0.1);%don't need to normalize by such a small amount /(1-normcdf(-11,0,1))/(1-normcdf(-123,0,1))

end