function kernel = computeFilter(theta,N)
    % returns convolution kernel function handle
    Kq = build_Kq([theta(1),theta(1)],N);
    AN = build_AN(Kq,theta(2),N);
    [ANhat,dANhat_dqM] = build_expm_stuff(AN,tau,1);
    BN=build_BN(theta(2),N);
    BNhat=build_BNhat(AN,ANhat,BN);
    kernel = @(s) ANhat^(floor(s/tau)-1)*BNhat;
end