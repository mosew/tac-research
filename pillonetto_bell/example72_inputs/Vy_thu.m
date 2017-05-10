function Vy_thu = Vy_thu(theta,n,u)
    
<<<<<<< HEAD
    
    tau = 1/n;
    Vy_thu = zeros(n);
    
=======
    tau = 1/n;
    conv_u = zeros(1,n);
        
>>>>>>> discrete_conv
    for i = 1:n
        Vy_thu(i,i) = (.05*L_i(theta,u,i,tau))^2;
    end
<<<<<<< HEAD
    
%     Vy_thu = .05*.05* diag(L_i_row(theta,u,tau)).^2;
    
    
%     Vy_thu = .05*eye(n);

%     for i = 1:n
%         for k = 1:i
%             Vy_thu(i,k) = conv_u(i)*conv_u(k);
%             if i~=k
%                 Vy_thu(k,i) = Vy_thu(i,k);
%             end
% 
%         end
%     end
%     Vy_thu = .05*.05*Vy_thu;
=======
        
    Vy_thu = .05*.05* diag(conv_u).^2;
    
>>>>>>> discrete_conv
end