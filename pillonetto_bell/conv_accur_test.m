theta = [exp(1),2];
tau = 1/50;

g = @(s) exp(-theta(2)*s);
f = @(s) 2*sin(s*8);

N = 4*10^3;
t = tau:(tau/N):1;

gsamp = exp(-theta(2).*t);
fsamp = feval(f,t);
output = conv(fsamp,gsamp,'full')*tau/N;
output = output(1:N:(50*N+1));

plot(output);
hold on
outty = zeros(1,50);
for i = 1:50
    fxg = @(s) feval(f,s).*feval(g,i*tau-s);
    outty(i) = integral(fxg,0,i*tau);
end
plot(outty)
%legend('Discrete','Continuous')