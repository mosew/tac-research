function output = L_i(theta,f,i,tau)
    % INPUTS:
    % i is the desired index of the output
    % theta is a vector of parameters
    % f is the input function
    % tau is the timestep
    
    % OUTPUTS:
    % output is a single number: it's the i^th output our system applied to f
    
    % FOR EXAMPLE 7.2 OF PILLONETTO-BELL
    
    g = @(s) exp(-theta(2)*s);
    try
        fxg = @(s) feval(f,s).*feval(g,i*tau-s);
    catch
        fxg = @(s) f.*feval(g,i*tau-s);
    end
    output = integral(fxg,0,i*tau);
end