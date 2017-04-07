function p_y_given_th = p_y_given_theta(y,th,tau,T,P,n,rkhs_eigenfile,data_path)
     p_y_given_th = mvnpdf(y,0,Vy_th(th,P,T,n,tau,rkhs_eigenfile,data_path));
end