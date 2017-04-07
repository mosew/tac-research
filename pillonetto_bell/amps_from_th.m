function a=amps_from_th(P,K,thetas,y,rkhs_eigenfile,T,n,tau,data_path)
    
    a = zeros(P,K);
    for k = 1:K
        EV=EV_aP_given_theta_y(thetas(:,k),y,rkhs_eigenfile,P,T,n,tau,data_path);
        mus = EV{1};
        covs = EV{2};
        a(:,k) = mvnrnd(mus,covs);
    end
    
end