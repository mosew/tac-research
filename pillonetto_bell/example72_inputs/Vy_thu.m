function Vy_thu = Vy_thu(theta,n,u)
    % I think this is where our covariance estimate will come in, but for
    % now we'll just use the stupid identity matrix. so STUPID
    Vy_thu = diag(.05*u);
end