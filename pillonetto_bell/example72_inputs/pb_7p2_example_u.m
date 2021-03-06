function input = pb_7p2_example_u()
    input = @(s) 0.3*betapdf(s,12,7) + 0.6*betapdf(s,4,11);
    input = @(s) feval(input,s) + normrnd(0,.05*feval(input,s));
end