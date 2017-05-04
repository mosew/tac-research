theta = [exp(1),10];
i = 15;
tau = 1/49;

g = @(s) exp(-theta(2)*s);
f = @(s) 2*sin(s);
fxg = @(s) feval(f,s).*feval(g,i*tau-s);

g_discrete = feval(g,tau*(0:(i-1)));
f_discrete = feval(f,tau*(0:(i-1)));
output = conv(f_discrete,g_discrete,'same');

[output(i), integral(fxg,0,i*tau,'RelTol',1e-3)]