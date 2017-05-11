function output = L_i(N,theta,f,i,tau)
    % INPUTS:
    % N is the spatial discretization
    % i is the desired index of the output
    % theta is a vector of parameters
    % f is the input function handle
    % tau is the timestep
    % dL_i_desired is a boolean which is 1 if our algorithm should
    %   concurrently compute dLmatrix for theta(1), which is efficiently done here
    %
    % OUTPUTS:
    % output is a single number: it's the i^th output of our system applied to f
    %
    % FOR ROSEN TRANSDERMAL ALCOHOL MODEL
    
    if ~isnumeric(f)        
        fxg = @(s) computeKernel(i*tau-s).*feval(f,s);
        output = integral(fxg,0,i*tau,'ArrayValued',true);
        output = output(1);        
    else
        % if f is a numeric vector, it's encoding amplitudes of kernel eigenfunctions
        fprintf('Somethin wrong with L_i\n');
    end

    
end