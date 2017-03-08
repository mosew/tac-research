function lo_up_mid = confidence_limits(fk)
    % fk is the input function after sampling a for the kth time.
    fL = @(t) quantile(feval(fk,t),0.025);
    fU = @(t) quantile(feval(fk,t),0.975);
    fM = @(t) mean(feval(fk,t));
    lo_up_mid={fL,fU,fM};
end