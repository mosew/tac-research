function A = acceptance(thetatry,thetaprev,y,tau,T,P,n,rkhs_eigenfile,data_path)
    % Computes acceptance ratio of draw
    % INPUT:
    % thetatry and thetaprev are nTheta x 1 vectors
    % 
    % OUTPUT:
    % single number
    
    assert(length(thetatry)==length(thetaprev));
    
<<<<<<< HEAD
    nTheta= length(thetatry);
    A=zeros(nTheta,1);
    num = p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path);
    den= p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path);
    for i = 1:nTheta
        p_thetatry = p_theta(i,thetatry);
        p_thetaprev = p_theta(i,thetaprev);
        A(i) = min(1, num*p_thetatry/(den*p_thetaprev));
        %A(i) = exp( log(p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path))+log(p_thetatry) - log(p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path))-log(p_thetaprev));
    end
=======
    p_thetatry=p_theta(thetatry);
    p_thetaprev=p_theta(thetaprev);
    num = p_y_given_theta(y,thetatry,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetatry;
    den = p_y_given_theta(y,thetaprev,tau,T,P,n,rkhs_eigenfile,data_path)*p_thetaprev;
    A = min(1, num/den);

>>>>>>> discrete_conv
end