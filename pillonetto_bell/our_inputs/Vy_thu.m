function Vy_thu = Vy_thu(N,theta,n,u)
    
    tau = 1/n;
    conv_u = zeros(1,n);
        
    for i = 1:n
        conv_u(i) = (.05*L_i(N,theta,u,i,tau)).^2;
    end
    
    Vy_thu = .05*.05* diag(conv_u).^2;
    
end