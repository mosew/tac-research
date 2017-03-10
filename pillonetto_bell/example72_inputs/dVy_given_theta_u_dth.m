function dVy_given_theta_u_dth = dVy_given_theta_u_dth(n)
    % This is the derivative of the known covariance matrix
    % V(y|theta,u)
    % with respect to theta
    % FOR EXAMPLE 7.2 FROM PILLONETTO-BELL
    
    % The paper is not totally clear to me, but I think the v_i are iid Gaussian
    % with mean 0 and variance .05u(t) at time t
    
    % In other words, V(y|theta,u) does not depend on theta.

    dVy_given_theta_u_dth = cell(1,2);
    dVy_given_theta_u_dth{1} = zeros(n);
    dVy_given_theta_u_dth{2} = zeros(n);
    
end