function Vy_thu = Vy_thu(theta,n,u)
    
    tau = 1/n;
    conv_u = zeros(1,n);
        
    for i = 1:n
        conv_u(i) = L_i(theta,u,i,tau);
    end
        
    Vy_thu = .05*.05* diag(conv_u).^2;
    
end