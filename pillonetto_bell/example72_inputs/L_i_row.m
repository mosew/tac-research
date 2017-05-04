function output = L_i_row(theta,f,tau)
    % INPUTS:
    % theta is a vector of parameters
    % f is the input function handle
    % tau is the timestep
    %
    % OUTPUTS:
    % output is a single number: it's the i^th output of our system applied to f
    %
    % FOR EXAMPLE 7.2 OF PILLONETTO-BELL
  
    N = 10^2;
    t = 0:(tau/(N-1)):1;
    
    gsamp = exp(-theta(2).*t);
    fsamp = feval(f,t);
    output = conv(fsamp,gsamp,'same');
    output = output(1:(N-1):end);
%     
%     fxg = @(s) feval(f,s).*feval(g,i*tau-s);
%     output = integral(fxg,0,i*tau,'RelTol',1e-3);

%     g_discrete = feval(g,tau*(0:(i-1)));
%     f_discrete = feval(f,tau*(0:(i-1)));
%     output = sum(conv(f_discrete,g_discrete,'same'));
    
end