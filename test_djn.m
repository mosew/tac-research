JN_and_dJN_globals
global M n
p0=ones(1,M+2+n);
[jn0,djn0]=JN_and_dJN(p0);

epsilon=10^(-6);
delta_jn = zeros(1,M+2+n);

for k=1:(M+2+n)
    p1 = zeros(1,M+2+n);
    p1(k) = epsilon;
    delta_jn(k)=(JN_and_dJN(p0+p1)-jn0)/epsilon;
end

% MSE componentwise
delta_jn-djn0
norm(delta_jn-djn0)