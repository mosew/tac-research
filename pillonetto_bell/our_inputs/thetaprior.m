function prior = thetaprior(theta)
    prior = normpdf(theta,.0046, 6e-4);
end