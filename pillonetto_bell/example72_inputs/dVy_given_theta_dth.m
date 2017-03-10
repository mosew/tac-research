function dVy_given_theta_dth = dVy_given_theta_dth(theta,y)
    % This one's gonna be for example 7.2 from the paper
    % Should be n x n, where n is the length of y
    
    n = length(y);
    
    dVy_given_theta_dth = cell(1,2);
    dVy_given_theta_u_dth = dVy_given_theta_u_dth(n);
    
    dVy_given_theta_dth{1} = Vy_given_theta/theta(1) - Vy_given_theta_u + dVy_given_theta_u_dth{1};
    dVy_given_theta_dth{2} = -2*theta(2)* Vy_given_theta - Vy_given_theta_u + dVy_given_theta_u_dth{2};
    
end