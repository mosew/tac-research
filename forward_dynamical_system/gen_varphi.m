function varphi = gen_varphi(ANhat,BNhat,u)
    global N n
    varphi=zeros(N+1,n+1);
    % goes from t = 0 to t = tau*n, where t = (j-1)*tau --> j = floor(t/tau) + 1
    for j=2:(n+1)
        varphi(:,j)=ANhat * varphi(:,j-1) + BNhat * u(j-1);
    end
end