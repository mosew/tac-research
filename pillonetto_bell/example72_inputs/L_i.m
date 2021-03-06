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
    
    if ~isnumeric(f)
    
        g = @(s) exp(-theta(2)*s);
        fxg = @(s) feval(f,s).*feval(g,i*tau-s);
        output = integral(fxg,0,i*tau);
        
    else
        % if f is a numeric vector, it's encoding amplitudes of kernel eigenfunctions
        % In the below, it's assumed that T=1
        P = length(f);
        w = pi/2*(2*(1:P)-1);
        output = sum( f.*sqrt(2)./(theta(2)^2+w.^2) .* (theta(2)*sin(w*i*tau) - w.*(cos(w*i*tau)) + w.*(exp(-theta(2)*i*tau))));
    end

    
end