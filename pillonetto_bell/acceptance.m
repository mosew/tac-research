function A = acceptance(thetatry,thetaprev,y,tau,T,P,n,rkhs_eigenfile,data_path)
    % Computes acceptance ratio of draw
    % for theta_i, the i^th parameter
    % INPUT:
    % thetatry and thetaprev are nTheta x 1 vectors
    % 
    % OUTPUT:
    % nTheta x 1 vector of numbers
    
    assert(length(thetatry)==length(thetaprev));
    
    nTheta= length(thetatry);
    A=zeros(nTheta,1);
    for i = 1:nTheta
        p_thetatry = p_theta(i,thetatry);
        p_thetaprev = p_theta(i,thetaprev);
        %A(i) = min(1, p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetatry / (p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetaprev));
        A(i) = exp( log(p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path))+log(p_thetatry) - log(p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path))-log(p_thetaprev));

    end
end