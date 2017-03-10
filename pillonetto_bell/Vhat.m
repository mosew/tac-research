function Vhat = Vhat(a)
    % This should compute the Cramer-Rao lower bound for our estimates of
    % EACH INDIVIDUAL PARAMETER
    %
    % a is the number of parameters
    % y is the output signal
    
    Vyth=Vy_given_theta;
    
    info = zeros(a);
    
    for i=1:a
        for j=1:i
            info(i,j) = 1/2 *trace(Vyth\dVy_given_theta_dth(i) * Vyth\dVy_given_theta_dth(j));
            if i~=j
                info(j,i) = info(i,j);
            end
        end
    end    
    
    Vhat = inv(info);
    
end