function output = L_i(N,theta,f,i,tau)
    % INPUTS:
    % i is the desired index of the output
    % theta is a vector of parameters
    % f is the input function handle
    % tau is the timestep
    %
    % OUTPUTS:
    % output is a single number: it's the i^th output of our system applied to f
    %
    % FOR ROSEN TRANSDERMAL ALCOHOL MODEL
    
    if ~isnumeric(f)
        % A little weird since we are doing M=0.
        
        g = computeFilter(theta,N);
        fxg = @(s) feval(f,s).*feval(g,i*tau-s);
        output = integral(fxg,0,i*tau);
        
    else
        % if f is a numeric vector, it's encoding amplitudes of kernel eigenfunctions
        fprintf('Somethin wrong with L_i');
    end

    
end