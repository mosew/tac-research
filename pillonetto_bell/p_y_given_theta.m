function p_y_given_th = p_y_given_theta(y,th,tau,T,P,rkhs_eigenfile,data_path)
    p_y_given_th = normpdf(y,0,Vy_given_theta(th,tau,T,P,rkhs_eigenfile,data_path));
%    p_y_given_th = det(2*pi*Vy_given_th(th,P))^-.5 * exp(-0.5 * y' * (Vy_given_th(th,P)\y));
end