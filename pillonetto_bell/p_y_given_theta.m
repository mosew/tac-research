function p_y_given_th = p_y_given_theta(y,th,tau,T,P,n,rkhs_eigenfile,data_path)
     p_y_given_th = mvnpdf(y,0,Vy_th(th,P,T,n,tau,rkhs_eigenfile,data_path));
%    p_y_given_th = det(2*pi*Vy_given_th(th,P))^-.5 * exp(-0.5 * y' * (Vy_given_th(th,P)\y));
end