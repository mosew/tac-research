function [eiv,eif] = haarwavelet_eigen(j,T)
    lmax = 0;
    lmin = -log2(T);
    d = str2double(dec2bin(j));
    l = -lmin - find(d);
    haar = @(t) (t>=0).*(t<1/2) - (t>=1/2).*(t<1);
    eiv = @(th) th(3)./c^(2*j);
    eif = @(s) c^(j/2).*haar(c^j .* s - k);
end