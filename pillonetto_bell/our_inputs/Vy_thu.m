function Vy_thu = Vy_thu(theta,n,u_total)
    % GOTTA FIGURE THIS OUT FOR ROSEN TRANSDERMAL ALCOHOL MODEL
    
    Vy_thu = cov(u_total);
    
%     tau = 1/n;
%     conv_u = zeros(1,n);
%         
%     for i = 1:n
%         conv_u(i) = L_i(N,theta,u,i,tau);
%     end
%     
%     Vy_thu = .1*.1* diag(conv_u.^2);
%     
end