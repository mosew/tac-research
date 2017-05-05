function Vy_thu = Vy_thu(theta,n,u)
    % I think this is where our covariance estimate will come in
    
    
    tau = 1/(n-1);
    Vy_thu = zeros(n);
    
    for i = 0:(n-1)
        Vy_thu(i+1,i+1) = (.05*L_i(theta,u,i,tau))^2;
    end
    
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
end