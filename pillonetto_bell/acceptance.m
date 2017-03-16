function A = acceptance(i,thetatry,thetaprev,y)
    % Computes acceptance ratio of draw
    % for theta_i, the i^th parameter
    p_thetatry = thetaprior(i,thetatry);
    p_thetaprev = thetaprior(i,thetaprev);
    A = min(1, p_y_given_theta(y,thetatry)*p_thetatry / (p_y_given_theta(y,thetaprev)*p_thetaprev));
end