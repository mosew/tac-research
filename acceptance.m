function A = acceptance(thetatry,thetaprev,y,tau,T,P,n,rkhs_eigenfile,data_path)
    % Computes acceptance ratio of draw
    % INPUT:
    % thetatry and thetaprev are nTheta x 1 vectors
    % 
    % OUTPUT:
    % single number
    
    assert(length(thetatry)==length(thetaprev));
    
    p_thetatry=p_theta(thetatry);
    p_thetaprev=p_theta(thetaprev);
    num = p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetatry;
    den = p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetaprev;
    A = min(1, num/den);

end