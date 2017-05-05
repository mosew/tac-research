theta = [exp(1),2];
i = 15;
tau = 1/49;

g = @(s) exp(-theta(2)*s);
f = @(s) 2*sin(s);
fxg = @(s) feval(f,s).*feval(g,i*tau-s);

N = 10^3;
t = 0:(tau/(N-1)):1;

gsamp = exp(-theta(2).*t);
fsamp = feval(f,t);
output = conv(fsamp,gsamp,'same');
output = output(1:(N-1):end);

plot(output);
hold on
outty = zeros(1,50);
for i = 1:50
    outty(i) = integral(fxg,0,(i-1)*tau);
end
plot(outty)
legend('Discrete','Continuous')