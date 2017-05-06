function Vy_thu = Vy_thu(theta,n,u)
    % I think this is where our covariance estimate will come in
    
    %Vy_thu = .05*eye(n);
    
    tau = 1/n;
    conv_u = zeros(1,n);
        
    for i = 1:n
        conv_u(i) = L_i(theta,u,i,tau);
    end
    
    %Vy_thu = .05 *.05 * (conv_u' * conv_u);
    
    Vy_thu = .05*.05* diag(conv_u).^2;
    
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
end