function K2 = greenKernel2(T)
    K2 = @(th,s,t) th*(s*t*min(s,t) - (s+t)*min(s,t).^2/2 + (min(s,t)).^3);
end