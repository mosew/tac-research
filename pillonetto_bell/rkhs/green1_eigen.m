function [eiv,eif] = green1_eigen(j,T)
    eiv = @(th) th(3)*T^2 / ( (j-1)*pi + pi/2)^2;
    eif = @(s) sqrt(2/T)*sin((s/T)*(j*pi-pi/2));
end