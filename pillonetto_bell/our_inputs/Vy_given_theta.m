function Vy_given_theta = Vy_given_theta(theta,T,P)
    % INPUTS:
    % theta is our parameter vector
    % T is the total length of an episode, in minutes
    % P is number of eigenfunctions to use
    
    data_path = '022717_234_splhr_arrays.mat';
    
    
    global tau n
    load(data_path);
    M=1;
    T=n*tau;
    q2=theta(1);
    qM=theta(2:end);
    t = 0:tau:T;
    
    % RKHS dependence
    [eivs,eifs] = get_kernel_eigenstuff(theta,P,T,'green2_eigen');
    sampled_eifs = sample_eigenfunctions(eifs,t,P);
    
    % Data dependence
    meas_noise = Vy_given_qu(data_path);
    
    
    convolved_sampled_eifs = zeros(size(sampled_eifs));
    for i = 1:P
        convolved_sampled_eifs(i,:) = CNhat*forward_system(q2,qM,tau,M,N,sampled_eifs(i,:));
    end
    
    
    Vy_given_theta = zeros(n+1,n+1);
    for i=1:(n+1)
        for k=1:(n+1)
            Vy_given_theta(i,k) = sum(eivs'.*convolved_sampled_eifs(:,i).*convolved_sampled_eifs(:,k)) + meas_noise(i,k);
        end
    end
    
    % Doesn't come out right. Maybe the episodes are too long?
    
end