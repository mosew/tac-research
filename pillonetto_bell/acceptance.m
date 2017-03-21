function A = acceptance(thetatry,thetaprev,y,tau,T,P,n,rkhs_eigenfile,data_path)
    % Computes acceptance ratio of draw
    % for theta_i, the i^th parameter
    assert(length(thetatry)==length(thetaprev));
    
    for i = 1:length(thetatry)
        p_thetatry = p_theta(i,thetatry);
        p_thetaprev = p_theta(i,thetaprev);
        A = min(1, p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetatry / (p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetaprev));
    end
end