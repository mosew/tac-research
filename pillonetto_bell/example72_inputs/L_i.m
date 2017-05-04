function output = L_i(theta,f,i,tau)
    % INPUTS:
    % i is the desired index of the output
    % theta is a vector of parameters
    % f is the input function handle
    % tau is the timestep
    %
    % OUTPUTS:
    % output is a single number: it's the i^th output of our system applied to f
    %
    % FOR EXAMPLE 7.2 OF PILLONETTO-BELL
    
    g = @(s) exp(-theta(2)*s);
    fxg = @(s) feval(f,s).*feval(g,i*tau-s);
    output = integral(fxg,0,i*tau,'RelTol',1e-3);

%     g_discrete = feval(g,tau*(0:(i-1)));
%     f_discrete = feval(f,tau*(0:(i-1)));
%     output = sum(conv(f_discrete,g_discrete,'same'));
    
end