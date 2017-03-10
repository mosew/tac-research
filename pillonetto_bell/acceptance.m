function A = acceptance(thetatry,thetaprev,y)
    % Computes acceptance ratio of draw
    A = min(1, p_y_given_theta(y,thetatry)*thetaprior(thetatry) / (p_y_given_theta(y,thetaprev)*thetaprior(thetaprev)));
end