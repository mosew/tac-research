function [eiv,eif] = green2_eigen(j,T)
    % This is for the second-order Green kernel
    % See appendix of paper for k=2 numerical scheme
    a=alpha(j);
    c3 = @(c4) c4*(2/(1+exp(-2*a))-1)/sin(a);
    c2 = @(c4) c4 - c3(c4)*exp(-a);
    c1 = @(c4) -c4 - c3(c4)*exp(-a);
    c4 = integral(@(t) (c1(1)*cos(a*t/T) + c2(1)*sin(a*t/T) + c3(1)*exp(-a*(T-t)/T) + exp(-a*t/T)).^2,0,T).^(-.5);
    
    eif = @(t) c1(c4) * cos(a*t/T)+ c2(c4)*sin(a*t/T) + c3(c4)*exp(-a*(T-t)/T) + c4*exp(-a*t/T);
    eiv = @(th) th(1)*(T/a)^4;
end