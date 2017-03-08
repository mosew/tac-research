% execute MCMC


% max number of iterations
Kmax = 200;

% initialize empty arrays to hold information about each iteration
q2_ks = zeros(1,Kmax);
q1_ks = zeros(1,Kmax); %constant k
theta = zeros(1,Kmax); %reproducing kernel parameter

% initialize parameter values
q2_ks(1) = 