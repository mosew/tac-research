function Vy_thu = Vy_thu(theta,n,u)
    % I think this is where our covariance estimate will come in
    
%     Vy_thu = diag(theta(2)^2*ones(1,n));
    Vy_thu = .01*eye(n);

end